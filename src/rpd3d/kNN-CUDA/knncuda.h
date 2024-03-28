// ninwang
bool power_dist_cuda_global_dev(const float* ref_dev, const int ref_nb,
                                const size_t ref_pitch, const float* ref_weight,
                                const float* query_dev, const int query_nb,
                                const size_t query_pitch, const int dim,
                                float* pdist_dev,
                                const size_t pdist_pitch_in_bytes);

// ninwang
bool knn_weighted_cuda_global_dev(const float* ref_dev, const int ref_nb,
                                  const size_t ref_pitch,
                                  const float* ref_weight,
                                  const float* query_dev, const int query_nb,
                                  const size_t query_pitch, const int dim,
                                  const int k, int* knn_index_dev,
                                  const size_t index_pitch);

// no use
bool knn_cuda_global_dev(const float* ref_dev, const int ref_nb,
                         const size_t ref_pitch, const float* query_dev,
                         const int query_nb, const size_t query_pitch,
                         const int dim, const int k, int* knn_index_dev,
                         const size_t index_pitch);

/**
 * For each input query point, locates the k-NN (indexes and distances) among
 * the reference points. This implementation uses global memory to store
 * reference and query points.
 *
 * @param ref        refence points
 * @param ref_nb     number of reference points
 * @param query      query points
 * @param query_nb   number of query points
 * @param dim        dimension of points
 * @param k          number of neighbors to consider
 * @param knn_dist   output array containing the query_nb x k distances
 * @param knn_index  output array containing the query_nb x k indexes
 */
bool knn_cuda_global(const float* ref, const int ref_nb, const float* query,
                     const int query_nb, const int dim, const int k,
                     float* knn_dist, int* knn_index);

/**
 * For each input query point, locates the k-NN (indexes and distances) among
 * the reference points. This implementation uses texture memory for storing
 * reference points  and memory to store query points.
 *
 * @param ref        refence points
 * @param ref_nb     number of reference points
 * @param query      query points
 * @param query_nb   number of query points
 * @param dim        dimension of points
 * @param k          number of neighbors to consider
 * @param knn_dist   output array containing the query_nb x k distances
 * @param knn_index  output array containing the query_nb x k indexes
 */
bool knn_cuda_texture(const float* ref, const int ref_nb, const float* query,
                      const int query_nb, const int dim, const int k,
                      float* knn_dist, int* knn_index);

// [no use]
// /**
//  * For each input query point, locates the k-NN (indexes and distances) among
//  * the reference points. Using cuBLAS, the computation of the distance matrix
//  * can be faster in some cases than other implementations despite being more
//  * complex.
//  *
//  * @param ref        refence points
//  * @param ref_nb     number of reference points
//  * @param query      query points
//  * @param query_nb   number of query points
//  * @param dim        dimension of points
//  * @param k          number of neighbors to consider
//  * @param knn_dist   output array containing the query_nb x k distances
//  * @param knn_index  output array containing the query_nb x k indexes
//  */
// bool knn_cublas(const float* ref, const int ref_nb, const float* query,
//                 const int query_nb, const int dim, const int k, float*
//                 knn_dist, int* knn_index);
