#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <OpenCL/cl.h>
#include <omp.h>
#include <utility>  // para std::move

#define BITS_PER_WORD 64  // Cada palabra de 64 bits (uint64_t)

using namespace std;

// Verifica errores de OpenCL
void check_error(cl_int err, const char* operation) {
    if (err != CL_SUCCESS) {
        cerr << "Error durante '" << operation << "': " << err << endl;
        exit(EXIT_FAILURE);
    }
}

// Carga los nombres de todas las muestras y cuenta kmers
void load_sample_names_and_num_kmers(
    const char* filename,
    vector<string>& sample_names,
    int& num_kmers
) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir el archivo '" << filename << "'." << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    // Leer cabecera
    if (getline(file, line)) {
        stringstream ss(line);
        string token;
        // Omitir primera columna (ID de kmer)
        getline(ss, token, '\t');
        // Leer nombres de muestras
        while (getline(ss, token, '\t')) {
            token.erase(remove(token.begin(), token.end(), '\r'), token.end());
            token.erase(remove(token.begin(), token.end(), '\n'), token.end());
            sample_names.push_back(token);
        }
    } else {
        cerr << "Error al leer la cabecera en '" << filename << "'." << endl;
        exit(EXIT_FAILURE);
    }
    // Contar kmers
    num_kmers = 0;
    while (getline(file, line)) {
        ++num_kmers;
    }
    file.close();
}

// Filtra muestras vacías: actualiza sample_names a solo las válidas
// y llena valid_indices con los índices originales de las muestras no vacías
void filter_empty_samples(
    const char* filename,
    vector<string>& sample_names,
    vector<int>& valid_indices
) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir '" << filename << "' para filtrar muestras vacías." << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    // Leer cabecera
    vector<string> header_names;
    if (getline(file, line)) {
        stringstream ss(line);
        string token;
        getline(ss, token, '\t');  // Omitir ID de kmer
        while (getline(ss, token, '\t')) {
            token.erase(remove(token.begin(), token.end(), '\r'), token.end());
            token.erase(remove(token.begin(), token.end(), '\n'), token.end());
            header_names.push_back(token);
        }
    } else {
        cerr << "Error al leer cabecera para filtrar muestras vacías." << endl;
        exit(EXIT_FAILURE);
    }
    int original_n = header_names.size();
    vector<bool> non_empty(original_n, false);
    // Recorrer cada fila de datos
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        getline(ss, token, '\t');  // Omitir ID de kmer
        for (int i = 0; i < original_n; ++i) {
            if (!getline(ss, token, '\t')) break;
            if (stoi(token) != 0) non_empty[i] = true;
        }
    }
    file.close();
    // Construir valid_indices y nuevo sample_names
    valid_indices.clear();
    vector<string> filtered_names;
    for (int i = 0; i < original_n; ++i) {
        if (non_empty[i]) {
            valid_indices.push_back(i);
            filtered_names.push_back(header_names[i]);
        }
    }
    sample_names.swap(filtered_names);
}

// Carga un bloque de datos (muestras originales) y empaqueta en bits
void load_data_block(
    const char* filename,
    const vector<int>& original_indices,
    int num_kmers,
    int num_words_per_sample,
    vector<uint64_t>& data_matrix_block
) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir '" << filename << "' en load_data_block." << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    // Leer cabecera y descartar
    getline(file, line);
    int block_n = original_indices.size();
    data_matrix_block.assign((size_t)block_n * num_words_per_sample, 0ULL);
    int kmer_idx = 0;
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        getline(ss, token, '\t');  // ID de kmer
        // Leer todos los valores de la fila
        vector<string> values;
        while (getline(ss, token, '\t')) values.push_back(token);
        int word_pos = kmer_idx / BITS_PER_WORD;
        int bit_off  = kmer_idx % BITS_PER_WORD;
        for (int s = 0; s < block_n; ++s) {
            int orig_idx = original_indices[s];
            int val = stoi(values[orig_idx]);
            if (val != 0) {
                data_matrix_block[(size_t)s * num_words_per_sample + word_pos] |= (1ULL << bit_off);
            }
        }
        ++kmer_idx;
    }
    file.close();
    if (kmer_idx != num_kmers) {
        cerr << "Discrepancia en número de kmers: leído=" << kmer_idx << " esperado=" << num_kmers << endl;
        exit(EXIT_FAILURE);
    }
}

// Kernel OpenCL para Jaccard (popcount)
const char* kernel_source = R"CLC(
__kernel void jaccard_kernel(
    __global const ulong* A,
    __global const ulong* B,
    __global float* D,
    const int num_words,
    const int nA,
    const int nB,
    const int sameBlock
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    if (i >= nA || j >= nB) return;
    if (sameBlock && i == j) {
        D[i * nB + j] = 0.0f;
        return;
    }
    int inter = 0, uni = 0;
    for (int w = 0; w < num_words; ++w) {
        ulong wa = A[i * num_words + w];
        ulong wb = B[j * num_words + w];
        ulong pi = wa & wb;
        ulong pu = wa | wb;
        inter += popcount(pi);
        uni   += popcount(pu);
    }
    float dist = (uni == 0) ? 1.0f : 1.0f - ((float)inter / (float)uni);
    D[i * nB + j] = dist;
}
)CLC";

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Uso: " << argv[0] << " <input> <output> [num_cpus]" << endl;
        return 1;
    }
    const char* input_file  = argv[1];
    const char* output_file = argv[2];
    int num_cpus = (argc >= 4) ? atoi(argv[3]) : omp_get_max_threads();
    omp_set_num_threads(num_cpus);

    time_t t0 = time(NULL);
    // 1) Cargar nombres y kmers
    vector<string> sample_names;
    int num_kmers;
    load_sample_names_and_num_kmers(input_file, sample_names, num_kmers);
    int original_n = sample_names.size();
    // 2) Filtrar muestras vacías
    vector<int> valid_indices;
    filter_empty_samples(input_file, sample_names, valid_indices);
    int n = sample_names.size();  // número de muestras tras filtrar
    cout << "Muestras originales: " << original_n
         << ", tras filtrar: " << n
         << ", kmers: " << num_kmers << endl;
    // 3) Construir mapeo original -> filtrado
    vector<int> to_filtered(original_n, -1);
    for (int f = 0; f < (int)valid_indices.size(); ++f) {
        to_filtered[ valid_indices[f] ] = f;
    }
    // 4) Precalcular palabras por muestra
    int num_words = (num_kmers + BITS_PER_WORD - 1) / BITS_PER_WORD;
    // 5) Preparar OpenCL
    cl_int err;
    cl_platform_id platform;
    err = clGetPlatformIDs(1, &platform, NULL);
    check_error(err, "clGetPlatformIDs");
    cl_device_id device;
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (err != CL_SUCCESS) {
        cout << "GPU no encontrada, usando CPU" << endl;
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL);
        check_error(err, "clGetDeviceIDs CPU");
    }
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    check_error(err, "clCreateContext");
    cl_program program = clCreateProgramWithSource(context, 1, &kernel_source, NULL, &err);
    check_error(err, "clCreateProgramWithSource");
    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t logsz;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsz);
        vector<char> logbuf(logsz);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logsz, logbuf.data(), NULL);
        cerr << "Build log:\n" << logbuf.data() << endl;
        return 1;
    }
    cl_kernel kernel = clCreateKernel(program, "jaccard_kernel", &err);
    check_error(err, "clCreateKernel");

    // 6) Dividir muestras en bloques secuenciales
    int batch = n; // sin partición: todo en uno (puedes ajustar)
    vector<vector<int>> blocks;
    for (int start = 0; start < n; start += batch) {
        int end = min(start + batch, n);
        vector<int> fb;
        for (int f = start; f < end; ++f) fb.push_back(f);
	blocks.push_back(std::move(fb));
    }
    int B = blocks.size();
    // 7) Matriz global de distancias inicializada en -1
    vector<float> G((size_t)n * n, -1.0f);

    // 8) Calcular pares de bloques (simétrico)
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int bi = 0; bi < B; ++bi) {
        for (int bj = bi; bj < B; ++bj) {
            auto &fbi = blocks[bi];
            auto &fbj = blocks[bj];
            int ni = fbi.size(), nj = fbj.size();
            // Construir índices originales para cada bloque
            vector<int> orig_i(ni), orig_j(nj);
            for (int ii = 0; ii < ni; ++ii) orig_i[ii] = valid_indices[fbi[ii]];
            for (int jj = 0; jj < nj; ++jj) orig_j[jj] = valid_indices[fbj[jj]];
            // Cargar datos
            vector<uint64_t> A, Bm;
            load_data_block(input_file, orig_i, num_kmers, num_words, A);
            load_data_block(input_file, orig_j, num_kmers, num_words, Bm);
            size_t szA = sizeof(uint64_t) * A.size();
            size_t szB = sizeof(uint64_t) * Bm.size();
            cl_mem bufA = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, szA, A.data(), &err);
            check_error(err, "clCreateBuffer A");
            cl_mem bufB = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, szB, Bm.data(), &err);
            check_error(err, "clCreateBuffer B");
            size_t pair_count = (size_t)ni * nj;
            vector<float> D(pair_count);
            cl_mem bufD = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * pair_count, NULL, &err);
            check_error(err, "clCreateBuffer D");
            // Argumentos
            err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufA); check_error(err, "arg0");
            err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufB); check_error(err, "arg1");
            err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &bufD); check_error(err, "arg2");
            err = clSetKernelArg(kernel, 3, sizeof(int), &num_words);       check_error(err, "arg3");
            err = clSetKernelArg(kernel, 4, sizeof(int), &ni);              check_error(err, "arg4");
            err = clSetKernelArg(kernel, 5, sizeof(int), &nj);              check_error(err, "arg5");
            int same = (bi==bj); err = clSetKernelArg(kernel, 6, sizeof(int), &same); check_error(err, "arg6");
            // Ejecutar
            cl_command_queue queue = clCreateCommandQueue(context, device, 0, &err);
            size_t gw[2] = {(size_t)ni, (size_t)nj};
            err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, gw, NULL, 0, NULL, NULL); check_error(err, "enqueue");
            clFinish(queue);
            err = clEnqueueReadBuffer(queue, bufD, CL_TRUE, 0, sizeof(float)*pair_count, D.data(), 0, NULL, NULL);
            check_error(err, "readD");
            // Volcar resultados a G usando índices filtrados
            for (int ii = 0; ii < ni; ++ii) {
                int gi = fbi[ii];
                for (int jj = 0; jj < nj; ++jj) {
                    int gj = fbj[jj];
                    float dist = D[(size_t)ii * nj + jj];
                    G[(size_t)gi * n + gj] = dist;
                    G[(size_t)gj * n + gi] = dist;
                }
            }
            // Liberar
            clReleaseMemObject(bufD);
            clReleaseMemObject(bufB);
            clReleaseMemObject(bufA);
            clReleaseCommandQueue(queue);
        }
    }
    // 9) Verificar pares no calculados
    int missing = 0;
    for (float d : G) if (d < 0.0f) ++missing;
    cout << "Pares con distancia no calculada: " << missing << " de " << (n*n) << endl;
    
    // 10) Escribir salida
    ofstream out(output_file);
    if (!out.is_open()) { cerr<<"Error abriendo salida."<<endl; return 1; }
    out << '\t';
    for (int i = 0; i < n; ++i) {
        out << sample_names[i] << (i+1<n? '\t':'\n');
    }
    for (int i = 0; i < n; ++i) {
        out << sample_names[i] << '\t';
        for (int j = 0; j < n; ++j) {
            out << G[(size_t)i*n + j] << (j+1<n? '\t':'\n');
        }
    }
    out.close();
    
    // Limpieza
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseContext(context);
    
    time_t t1 = time(NULL);
    cout << "Tiempo total: " << (t1 - t0) << " s" << endl;
    return 0;
}

