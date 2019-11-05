#include <cstdio>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include "include/stopwatch.h"
#include "include/tweets_spatial_analysis.h"

// template void build_adjacency_matrix<short>(string input_file, string output_file, long rows = 0);
// template void find_components<short>(string adj_file, string outlier_file = "");
// template vector<pair<int, vector<int>>> bfs<short>(Mat<short> &A);
// template void build_distance_matrix<short>(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
// template Mat<short> APD_recursive<short>(const Mat<short> &A, int rank, int n_procs);
// template Mat<short> APD<short>(const Mat<short> &A, int rank, int n_procs);
// template Mat<short> BPWM<short>(const Mat<short> &A, const Mat<short> &B, int node, int n_procs);
// template Mat<short> compute_successor_matrix<short>(const Mat<short> &A, const Mat<short> &D, int node, int n_procs);
// template Mat<short> compute_successor_matrix_saving_witness<short>(const Mat<short> &A, const Mat<short> &D, int node, int n_procs, int i, string dis_file);
// template void build_successor_matrix<short>(string adj_file, string dis_file, int node, int n_procs);
// template void APSP<short>(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
// template void test_APD_recursive<short>(string mat_file, int mode, int node, int n_procs);
// template void test_APD<short>(string mat_file, int mode, int node, int n_procs);
//
// template void build_adjacency_matrix<int>(string input_file, string output_file, long rows = 0);
// template void find_components<int>(string adj_file, string outlier_file = "");
// template vector<pair<int, vector<int>>> bfs<int>(Mat<int> &A);
// template void build_distance_matrix<int>(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
// template Mat<int> APD_recursive<int>(const Mat<int> &A, int rank, int n_procs);
// template Mat<int> APD<int>(const Mat<int> &A, int rank, int n_procs);
// template Mat<int> BPWM<int>(const Mat<int> &A, const Mat<int> &B, int node, int n_procs);
// template Mat<int> compute_successor_matrix<int>(const Mat<int> &A, const Mat<int> &D, int node, int n_procs);
// template Mat<int> compute_successor_matrix_saving_witness<int>(const Mat<int> &A, const Mat<int> &D, int node, int n_procs, int i, string dis_file);
// template void build_successor_matrix<int>(string adj_file, string dis_file, int node, int n_procs);
// template void APSP<int>(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
// template void test_APD_recursive<int>(string mat_file, int mode, int node, int n_procs);
// template void test_APD<int>(string mat_file, int mode, int node, int n_procs);
//
// template void build_adjacency_matrix<long>(string input_file, string output_file, long rows = 0);
// template void find_components<long>(string adj_file, string outlier_file = "");
// template vector<pair<int, vector<int>>> bfs<long>(Mat<long> &A);
// template void build_distance_matrix<long>(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
// template Mat<long> APD_recursive<long>(const Mat<long> &A, int rank, int n_procs);
// template Mat<long> APD<long>(const Mat<long> &A, int rank, int n_procs);
// template Mat<long> BPWM<long>(const Mat<long> &A, const Mat<long> &B, int node, int n_procs);
// template Mat<long> compute_successor_matrix<long>(const Mat<long> &A, const Mat<long> &D, int node, int n_procs);
// template Mat<long> compute_successor_matrix_saving_witness<long>(const Mat<long> &A, const Mat<long> &D, int node, int n_procs, int i, string dis_file);
// template void build_successor_matrix<long>(string adj_file, string dis_file, int node, int n_procs);
// template void APSP<long>(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
// template void test_APD_recursive<long>(string mat_file, int mode, int node, int n_procs);
// template void test_APD<long>(string mat_file, int mode, int node, int n_procs);
//
// template void build_adjacency_matrix<double>(string input_file, string output_file, long rows = 0);
// template void find_components<double>(string adj_file, string outlier_file = "");
// template vector<pair<int, vector<int>>> bfs<double>(Mat<double> &A);
// template void build_distance_matrix<double>(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
// template Mat<double> APD_recursive<double>(const Mat<double> &A, int rank, int n_procs);
// template Mat<double> APD<double>(const Mat<double> &A, int rank, int n_procs);
// template Mat<double> BPWM<double>(const Mat<double> &A, const Mat<double> &B, int node, int n_procs);
// template Mat<double> compute_successor_matrix<double>(const Mat<double> &A, const Mat<double> &D, int node, int n_procs);
// template Mat<double> compute_successor_matrix_saving_witness<double>(const Mat<double> &A, const Mat<double> &D, int node, int n_procs, int i, string dis_file);
// template void build_successor_matrix<double>(string adj_file, string dis_file, int node, int n_procs);
// template void APSP<double>(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
// template void test_APD_recursive<double>(string mat_file, int mode, int node, int n_procs);
// template void test_APD<double>(string mat_file, int mode, int node, int n_procs);

template void build_adjacency_matrix<float>(string input_file, string output_file, long rows = 0);
template void find_components<float>(string adj_file, string outlier_file = "");
template vector<pair<int, vector<int>>> bfs<float>(Mat<float> &A);
template void build_distance_matrix<float>(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
template Mat<float> APD_recursive<float>(const Mat<float> &A, int rank, int n_procs);
template Mat<float> APD<float>(const Mat<float> &A, int rank, int n_procs);
template Mat<float> BPWM<float>(const Mat<float> &A, const Mat<float> &B, int node, int n_procs);
template Mat<float> get_witness_matrix(const Mat<float> &A, const Mat<float> &D, int node, int n_procs, int i, int r, string dis_file);
template void build_witness_matrix<float>(string adj_file, string dis_file, int r, int node, int n_procs);
template Mat<float> compute_successor_matrix<float>(const Mat<float> &A, const Mat<float> &D, int node, int n_procs, int i, string dis_file);
template void build_successor_matrix<float>(string adj_file, string dis_file, int node, int n_procs);
template void APSP<float>(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
template void test_APD_recursive<float>(string mat_file, int mode, int node, int n_procs);
template void test_APD<float>(string mat_file, int mode, int node, int n_procs);

template <class T>
void build_adjacency_matrix(string ellipse_file, string adj_file, long rows) {
  shino::precise_stopwatch stopwatch;
  mat A;
  vector<int> remove_idx;
  vector<pair<unsigned int, long unsigned>> id;

  A.load(ellipse_file, csv_ascii);
  if (rows == 0) rows = A.n_rows;
  Mat<short> Adj(rows, rows, fill::zeros);
  cout << ellipse_file << ": " << rows << endl;
  long m = rows * rows / 2, mm = 0;
  for (int i = 0; i < rows; ++i) {
    if (A(i, 2) == 0 || A(i, 3) == 0) {
      remove_idx.push_back(i);
      // Adj.row(i) = -1 * ones<ivec>(rows).t();
      continue;
    }
    for (int j = 0; j < i; ++j) {
      Ellipse e1 = Ellipse(A(i, 0), A(i, 1), A(i, 2), A(i, 3), A(i, 4));
      Ellipse e2 = Ellipse(A(j, 0), A(j, 1), A(j, 2), A(j, 3), A(j, 4));
      if (e1.overlap(e2))
        Adj(i, j) = 1;
      else
        Adj(i, j) = 0;
      ++mm;
      if (mm % 1000000 == 0) {
        printf("%7.4f%% completed\n", (double)mm / m * 100);
      }
    }
    id.push_back({(unsigned int)i, (long unsigned)A(i, 5)});
  }
  Adj = Adj + Adj.t();

  // remove rows and cols that contain -1
  reverse(remove_idx.begin(), remove_idx.end());
  for (auto &i : remove_idx) {
    Adj.shed_row(i);
    Adj.shed_col(i);
  }
  rows = Adj.n_rows;
  Mat<long unsigned> id_mat(rows, 2, fill::zeros);
  for (int i = 0; i < rows; ++i) {
    id_mat(i, 1) = id[i].second;
  }
  // put the degree of node in the 1st column of id_mat
  id_mat.col(0) = conv_to<Mat<long unsigned>>::from(sum(Adj, 1));

  // output
  string adj_id_file(adj_file);
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  id_mat.save(adj_id_file, csv_ascii);
  Adj.save(adj_file, csv_ascii);

  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
}

template <class T>
vector<pair<int, vector<int>>> bfs(Mat<T> &A) {
  int rows = A.n_rows;
  Col<T> visited(rows, fill::zeros);
  vector<pair<int, vector<int>>> components;
  for (int i = 0; i < rows; ++i) {
    if (!visited(i)) {
      visited(i) = 1;
      vector<int> component({i});
      queue<int> q;
      for (int j = 0; j < rows; ++j) {
        if (!visited(j) && A(i, j) > 0) {
          visited(j) = 1;
          q.push(j);
          component.push_back(j);
        }
      }
      while (!q.empty()) {
        int j = q.front();
        q.pop();
        for (int k = 0; k < rows; ++k) {
          if (!visited(k) && A(j, k) > 0) {
            visited(k) = 1;
            component.push_back(k);
            q.push(k);
          }
        }
      }
      components.push_back({(int)component.size(), component});
    }
  }
  return components;
}

template <class T>
void find_components(string adj_file, string outlier_file) {
  shino::precise_stopwatch stopwatch;
  Mat<T> A;
  Mat<long unsigned> id_mat;
  A.load(adj_file, csv_ascii);
  string adj_id_file = adj_file;
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  id_mat.load(adj_id_file);

  // string adj_id_file_update = adj_id_file;
  // adj_id_file_update.replace(adj_id_file_update.end() - 4,
  // adj_id_file_update.end(), "_update.csv"); id_mat.col(0) = conv_to<Mat<long
  // unsigned>>::from(sum(A, 1)); id_mat.save(adj_id_file_update, csv_ascii);

  cout << outlier_file << endl;
  if (outlier_file.size() > 5) {
    Col<long unsigned> outlier_id;
    outlier_id.load(outlier_file);
    unordered_map<long unsigned, int> id_map;
    int n = A.n_rows;
    for (int i = 0; i < n; ++i) {
      id_map[id_mat(i, 1)] = i;
    }

    int m = outlier_id.n_rows;
    vector<int> outlier_idx;
    for (int i = 0; i < m; ++i) {
      outlier_idx.push_back(id_map[outlier_id(i)]);
    }

    sort(outlier_idx.begin(), outlier_idx.end());
    reverse(outlier_idx.begin(), outlier_idx.end());

    uvec indices(m);
    for (int i = 0; i < m; ++i) {
      indices(i) = outlier_idx[i];
    }

    A.shed_rows(indices);
    A.shed_cols(indices);
    id_mat.shed_rows(indices);

    adj_file.replace(adj_file.end() - 4, adj_file.end(), "_" + to_string(m) + "_outlier.csv");
  }

  vector<pair<int, vector<int>>> components = bfs(A);
  sort(components.begin(), components.end());
  reverse(components.begin(), components.end());

  int rows = A.n_rows;
  int cur = 0;
  vector<vector<long unsigned>> id;
  for (int i = 0; i < rows; i++) {
    id.push_back({(long unsigned)i, (long unsigned)id_mat(i, 1)});
  }
  Mat<short> B(rows, rows, fill::zeros);
  int k = 0;
  for (auto &c : components) {
    int m = c.first;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        B(cur + i, cur + j) = (short)A(c.second[i], c.second[j]);
      }
      id[c.second[i]][0] = cur + i;
      id[c.second[i]].insert(id[c.second[i]].begin() + 1, k);
    }
    cur += m;
    ++k;
  }
  sort(id.begin(), id.end());
  for (int i = 0; i < rows; ++i) {
    id_mat(i, 0) = id[i][1];
    id_mat(i, 1) = id[i][2];
  }

  // output
  adj_file.replace(adj_file.end() - 4, adj_file.end(), "_ordered.csv");
  adj_id_file = adj_file;
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  id_mat.save(adj_id_file, csv_ascii);
  B.save(adj_file, csv_ascii);

  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
}

template <class T>
Mat<T> APD_recursive(const Mat<T> &A, int node, int n_procs) {
  long n = A.n_rows;
  Mat<T> Z;

  if (n >= 100 && n_procs > 1) {
    Z = matrix_square(A, node, n_procs);
  } else {
    Z = A * A;
  }

  Mat<T> B(n, n, fill::zeros);
  long cnt = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j && (A(i, j) == 1 || Z(i, j) > 0)) {
        B(i, j) = 1;
        ++cnt;
      } else {
        B(i, j) = 0;
      }
    }
  }

  Mat<T> D(n, n, fill::ones);
  D = -1 * D;
  if (cnt == (n - 1) * n) {
    D = 2 * B - A;
    return D;
  }

  Mat<T> S = APD_recursive<T>(B, node, n_procs);

  Mat<T> X;
  if (n >= 100 && n_procs > 1) {
    X = matrix_multiplication(S, A, node, n_procs);
  } else {
    X = S * A;
  }

  vec deg(n, fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      deg(i) += A(i, j);
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (X(i, j) >= S(i, j) * deg(j)) {
        D(i, j) = 2 * S(i, j);
      } else {
        D(i, j) = 2 * S(i, j) - 1;
      }
    }
  }
  return D;
}

template <class T>
Mat<T> APD(const Mat<T> &AA, int node, int n_procs) {
  long n = AA.n_rows;

  stack<Mat<T>> A_stack;
  Mat<T> A = AA;
  Mat<T> D;
  while (true) {
    Mat<T> Z;

    if (n >= 100 && n_procs > 1) {
      Z = matrix_square(A, node, n_procs);
    } else {
      Z = A * A;
    }

    Mat<T> B(n, n, fill::zeros);
    long cnt = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j && (A(i, j) == 1 || Z(i, j) > 0)) {
          B(i, j) = 1;
          ++cnt;
        } else {
          B(i, j) = 0;
        }
      }
    }

    if (cnt == (n - 1) * n) {
      D = 2 * B - A;
      break;
    }

    A_stack.push(A);
    A = B;
  }
  while (!A_stack.empty()) {
    A = A_stack.top();
    A_stack.pop();

    Mat<T> S = D;
    Mat<T> X;
    if (n >= 100 && n_procs > 1) {
      X = matrix_multiplication(S, A, node, n_procs);
    } else {
      X = S * A;
    }

    vec deg(n, fill::zeros);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        deg(i) += A(i, j);
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (X(i, j) >= S(i, j) * deg(j)) {
          D(i, j) = 2 * S(i, j);
        } else {
          D(i, j) = 2 * S(i, j) - 1;
        }
      }
    }
  }
  return D;
}

template <class T>
void build_distance_matrix(string adj_file, string dis_file, int node, int n_procs, bool non_recursive) {
  shino::precise_stopwatch stopwatch;
  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();

  Mat<T> A;
  Mat<long unsigned> id_mat;
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_file.c_str(), elapsed_time);
  }
  A.load(adj_file, csv_ascii);

  int rows = A.n_rows;
  string adj_id_file = adj_file;
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_id_file.c_str(), elapsed_time);
  }
  id_mat.load(adj_id_file);

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Unpacking id... (%u ms)\n", node, elapsed_time);
  }
  vector<pair<unsigned int, long unsigned>> id;
  for (int i = 0; i < rows; i++) {
    id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
  }

  vector<pair<int, int>> idx;
  for (int i = 0; i < rows;) {
    int start = i;
    while (i < rows && id[start].first == id[i].first) {
      ++i;
    }
    int end = i - 1;
    idx.push_back({start, end});
  }

  Mat<int> AA = -1 * ones<Mat<int>>(rows, rows);
  for (int i = 0; i < (int)idx.size(); ++i) {
    int r1 = idx[i].first;
    int c1 = idx[i].first;
    int r2 = idx[i].second;
    int c2 = idx[i].second;

    if (node == 0) {
      elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
      printf("CPU %d: Applying APD to block matrix %d (%d-by-%d) (%u ms)\n", node, i, r2 - r1 + 1, c2 - c1 + 1, elapsed_time);
    }

    Mat<T> D;
    if (non_recursive) {
      D = APD<T>(A.submat(r1, c1, r2, c2), node, n_procs);
    } else {
      D = APD_recursive<T>(A.submat(r1, c1, r2, c2), node, n_procs);
      dis_file.replace(dis_file.end() - 4, dis_file.end(), "_recursive.csv");
    }
    if (node == 0) {
      for (int j = r1; j < r2 + 1; ++j) {
        for (int k = c1; k < c2 + 1; ++k) {
          AA(j, k) = (int)D(j - r1, k - c1);
        }
      }
    }
  }

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Writing to %s... (%u ms)\n", node, dis_file.c_str(), elapsed_time);
    AA.save(dis_file, csv_ascii);

    cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
  }
}

template <class T>
Mat<T> DPWMD(const Mat<T> &A, const Mat<T> &B, int node, int n_procs) {
  Mat<T> C = matrix_multiplication(A, B, node, n_procs);
  Mat<T> W;
  std::set<pair<int, int>> L;
  long n = A.n_rows;
  for (long i = 0; i < n; ++i)
    for (long j = 0; j < n; ++j)
      if (C(i, j) > 0) L.insert(make_pair(i, j));
  while (!L.empty()) {
    Mat<T> R(n, n, fill::ones);
    long u = ceil(3 * log2(n) / log2(4 / 3));
    for (int k = 0; k < u; ++k) {
    }
  }

  return W;
}

template <class T>
Mat<T> BPWM(const Mat<T> &A, const Mat<T> &B, int node, int n_procs) {
  shino::precise_stopwatch stopwatch;
  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  Mat<T> W = -1 * matrix_multiplication(A, B, node, n_procs);
  long n = A.n_rows;
  long l = ceil(log2(n));
  long u = ceil(3.77 * log2(n));
  for (int k = 0; k < l; ++k) {
    long d = pow(2, k);
    for (int m = 0; m < u; ++m) {
      arma_rng::set_seed_random();
      Col<int> K = randi<Col<int>>(d, distr_param(0, n - 1));
      Mat<T> X(n, d, fill::zeros), Y(d, n, fill::zeros);
      for (int i = 0; i < d; ++i) {
        X.col(i) = (K(i) + 1) * A.col(K(i));
        Y.row(i) = B.row(K(i));
      }

      if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf(
            "CPU %d: Repeating randomized part... (k, d, m: %d, %ld, %d, %u "
            "ms)\n",
            node, k, d, m, elapsed_time);
      }

      Mat<T> C;
      if (k >= 11) {
        C = matrix_multiplication(X, Y, node, n_procs);
      } else {
        C = X * Y;
      }

      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          if (W(i, j) < 0 && C(i, j) > 0 && C(i, j) < n && A(i, (int)C(i, j) - 1) > 0 && B((int)C(i, j) - 1, j) > 0) {
            W(i, j) = C(i, j) - 1;
          }
        }
      }
    }
  }

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf(
        "CPU %d: Using brutal force to calculate the undecided entries in W... "
        "(%u ms)\n",
        node, elapsed_time);
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (W(i, j) < 0) {
        for (int k = 0; k < n; ++k) {
          if (A(i, k) > 0 && B(k, j) > 0) {
            W(i, j) = k;
            break;
          }
        }
      }
    }
  }
  return W;
}

inline bool exists(const string &name) {
  ifstream f(name.c_str());
  return f.good();
}

template <class T>
Mat<T> get_witness_matrix(const Mat<T> &A, const Mat<T> &D, int node, int n_procs, int i, int r, string dis_file) {
  string suffix = "witness_matrix_" + to_string(i) + "_" + to_string(r) + ".csv";
  string wit_file(dis_file);
  wit_file.replace(wit_file.end() - 19, wit_file.end(), suffix);

  long rows = D.n_rows;
  long cols = D.n_cols;
  Mat<T> W;

  if (exists(wit_file)) {
    if (node == 0) {
      printf("CPU %d: r: %d, Reading BPWM (%d, %d) from %s\n", node, r, i, r, wit_file.c_str());
    }
    W.load(wit_file, csv_ascii);
    return W;
  }

  Mat<T> Dr(rows, cols, fill::zeros);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (((int)D(i, j) + 1) % 3 == r) {
        Dr(i, j) = 1;
      }
    }
  }

  W = BPWM(A, Dr, node, n_procs);
  if (rows > 1000 || cols > 1000) {
    if (node == 0) {
      W.save(wit_file, csv_ascii);
      printf("CPU %d: r: %d, Saving BPWM (%d, %d) to %s\n", node, r, i, r, wit_file.c_str());
    }
  }
  return W;
}

template <class T>
Mat<T> compute_successor_matrix(const Mat<T> &A, const Mat<T> &D, int node, int n_procs, int i, string dis_file) {
  shino::precise_stopwatch stopwatch;
  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  long rows = D.n_rows;
  long cols = D.n_cols;

  vector<Mat<T>> Wr;
  Wr.reserve(3);
  for (int r = 0; r < 3; ++r) {
    Mat<T> W = get_witness_matrix(A, D, node, n_procs, i, r, dis_file);
    Wr.push_back(W);
  }

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf(
        "CPU %d: Constructing successor matrix from witness matrix... (%u "
        "ms)\n",
        node, elapsed_time);
  }
  Mat<T> S(rows, cols, fill::zeros);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      S(i, j) = Wr[(int)D(i, j) % 3](i, j);
    }
  }
  return S;
}

template <class T>
void build_witness_matrix(string adj_file, string dis_file, int r, int node, int n_procs) {
  shino::precise_stopwatch stopwatch;
  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  Mat<T> A, D;
  Mat<long unsigned> id_mat;
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_file.c_str(), elapsed_time);
  }
  A.load(adj_file, csv_ascii);

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, dis_file.c_str(), elapsed_time);
  }
  D.load(dis_file, csv_ascii);

  int rows = A.n_rows;
  string adj_id_file = adj_file;
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_id_file.c_str(), elapsed_time);
  }
  id_mat.load(adj_id_file);

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Unpacking id... (%u ms)\n", node, elapsed_time);
  }
  vector<pair<unsigned int, long unsigned>> id;
  for (int i = 0; i < rows; i++) {
    id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
  }

  vector<pair<int, int>> idx;
  for (int i = 0; i < rows;) {
    int start = i;
    while (i < rows && id[start].first == id[i].first) {
      ++i;
    }
    int end = i - 1;
    idx.push_back({start, end});
  }

  for (int i = 0; i < (int)idx.size(); ++i) {
    int r1 = idx[i].first;
    int c1 = idx[i].first;
    int r2 = idx[i].second;
    int c2 = idx[i].second;

    if (r2 - r1 + 1 >= 1000 || c2 - c1 + 1 >= 1000) {
      if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf(
            "CPU %d: Computing the witness matrix of block matrix %d - %d "
            "(%d-by-%d) (%u ms)\n",
            node, i, r, r2 - r1 + 1, c2 - c1 + 1, elapsed_time);
      }
      Mat<T> W = get_witness_matrix<T>(A.submat(r1, c1, r2, c2), D.submat(r1, c1, r2, c2), node, n_procs, i, r, dis_file);
    }
    if (node == 0) {
      elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
      cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
    }
  }
}

template <class T>
void build_successor_matrix(string adj_file, string dis_file, int node, int n_procs) {
  shino::precise_stopwatch stopwatch;
  auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
  Mat<T> A, D;
  Mat<long unsigned> id_mat;
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_file.c_str(), elapsed_time);
  }
  A.load(adj_file, csv_ascii);

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, dis_file.c_str(), elapsed_time);
  }
  D.load(dis_file, csv_ascii);

  int rows = A.n_rows;
  string adj_id_file = adj_file;
  adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_id_file.c_str(), elapsed_time);
  }
  id_mat.load(adj_id_file);

  if (node == 0) {
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Unpacking id... (%u ms)\n", node, elapsed_time);
  }
  vector<pair<unsigned int, long unsigned>> id;
  for (int i = 0; i < rows; i++) {
    id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
  }

  vector<pair<int, int>> idx;
  for (int i = 0; i < rows;) {
    int start = i;
    while (i < rows && id[start].first == id[i].first) {
      ++i;
    }
    int end = i - 1;
    idx.push_back({start, end});
  }

  Mat<int> SS = -1 * ones<Mat<int>>(rows, rows);
  for (int i = 0; i < (int)idx.size(); ++i) {
    int r1 = idx[i].first;
    int c1 = idx[i].first;
    int r2 = idx[i].second;
    int c2 = idx[i].second;

    if (node == 0) {
      elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
      printf(
          "CPU %d: Computing the successor matrix of block matrix %d "
          "(%d-by-%d) (%u ms)\n",
          node, i, r2 - r1 + 1, c2 - c1 + 1, elapsed_time);
    }

    Mat<T> S = compute_successor_matrix<T>(A.submat(r1, c1, r2, c2), D.submat(r1, c1, r2, c2), node, n_procs, i, dis_file);
    if (node == 0) {
      for (int j = r1; j < r2 + 1; ++j) {
        for (int k = c1; k < c2 + 1; ++k) {
          SS(j, k) = (int)S(j - r1, k - c1) + r1;
        }
      }
      elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
      printf(
          "CPU %d: Finish computing the successor matrix of block matrix %d "
          "(%d-by-%d) (%u ms)\n",
          node, i, r2 - r1 + 1, c2 - c1 + 1, elapsed_time);
    }
  }

  if (node == 0) {
    string suc_file = dis_file;
    suc_file.replace(suc_file.end() - 19, suc_file.end(), "successor_matrix.csv");
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    printf("CPU %d: Writing to %s... (%u ms)\n", node, suc_file.c_str(), elapsed_time);
    SS.save(suc_file, csv_ascii);

    cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
  }
}

template <class T>
void APSP(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs) {
  build_adjacency_matrix<T>(ellipse_file, adj_file, rows);
  find_components<T>(adj_file);
  build_distance_matrix<T>(adj_ordered_file, dis_file, node, n_procs);
  build_successor_matrix<T>(adj_ordered_file, dis_file, node, n_procs);
}

template <class T>
void test_APD(string mat_file, int mode, int node, int n_procs) {
  Mat<T> A;
  A.load(mat_file, csv_ascii);
  Mat<T> D = APD<T>(A, node, n_procs);
  if (mode == 2) {
    D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel_" + to_string(node) + ".csv"), csv_ascii);
  } else {
    if (node == 0) {
      D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel.csv"), csv_ascii);
    }
  }
}

template <class T>
void test_APD_recursive(string mat_file, int mode, int node, int n_procs) {
  Mat<T> A;
  A.load(mat_file, csv_ascii);
  Mat<T> D = APD_recursive<T>(A, node, n_procs);
  if (mode == 2) {
    D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel_" + to_string(node) + ".csv"), csv_ascii);
  } else {
    if (node == 0) {
      D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel.csv"), csv_ascii);
    }
  }
}

vector<int> get_shortest_path_by_index(Mat<int> &S, Mat<long unsigned> &id_mat, int from, int to) {
  vector<int> index_path{from};

  if (id_mat(from, 0) != id_mat(to, 0)) return index_path;

  for (; from != to; from = S(from, to)) index_path.push_back(S(from, to));

  return index_path;
}

vector<long unsigned> convert_index_path_to_id_path(Mat<long unsigned> &id_mat, vector<int> index_path) {
  int n = index_path.size();
  vector<long unsigned> id_path(n, 0);

  for (int i = 0; i < n; ++i) id_path[i] = id_mat(index_path[i], 1);

  return id_path;
}

pair<vector<int>, vector<long unsigned>> get_shortest_path_by_id(Mat<int> &S, Mat<long unsigned> &id_mat, long unsigned id_from, long unsigned id_to) {
  unordered_map<long unsigned, int> id_map;
  int n = id_mat.n_rows;

  for (int i = 0; i < n; ++i) id_map[id_mat(i, 1)] = i;

  int from = id_map[id_from];
  int to = id_map[id_to];
  vector<int> index_path = get_shortest_path_by_index(S, id_mat, from, to);
  vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);

  return make_pair(index_path, id_path);
}

vector<vector<int>> find_all_shortest_index_paths(Mat<int> &S, Mat<long unsigned> &id_mat) {
  vector<vector<int>> index_paths;
  int n = id_mat.n_rows;

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      vector<int> index_path = get_shortest_path_by_index(S, id_mat, i, j);
      if (index_path.size() <= 1) continue;
      index_paths.push_back(index_path);
    }
  }

  return index_paths;
}

void convert_paths_to_json(vector<vector<int>> &index_paths) {
  cout << "{";
  cout << "\"" << index_paths[0].front() << "\": ";
  cout << "{\"" << index_paths[0].back() << "\": [";
  for (int j = 0; j < (int)index_paths[0].size() - 1; ++j) cout << "\"" << index_paths[0][j] << "\", ";
  cout << "\"" << index_paths[0].back() << "\"]";
  for (int i = 1; i < (int)index_paths.size(); ++i) {
    if (index_paths[i].front() == index_paths[i - 1].front()) {
      cout << ",\n";
      cout << "  \"" << index_paths[i].back() << "\": [";
    } else {
      cout << "},\n";
      cout << " \"" << index_paths[i].front() << "\": {";
      cout << "\"" << index_paths[i].back() << "\": [";
    }
    for (int j = 0; j < (int)index_paths[i].size() - 1; ++j) cout << "\"" << index_paths[i][j] << "\", ";
    cout << "\"" << index_paths[i].back() << "\"]";
  }
  cout << "}}\n";
}

void convert_paths_to_json(vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  fout << "{";
  fout << "\"" << index_paths[0].front() << "\": ";
  fout << "{\"" << index_paths[0].back() << "\": [";
  for (int j = 0; j < (int)index_paths[0].size() - 1; ++j) fout << "\"" << index_paths[0][j] << "\", ";
  fout << "\"" << index_paths[0].back() << "\"]";
  for (int i = 1; i < (int)index_paths.size(); ++i) {
    if (index_paths[i].front() == index_paths[i - 1].front()) {
      fout << ",\n";
      fout << "  \"" << index_paths[i].back() << "\": [";
    } else {
      fout << "},\n";
      fout << " \"" << index_paths[i].front() << "\": {";
      fout << "\"" << index_paths[i].back() << "\": [";
    }
    for (int j = 0; j < (int)index_paths[i].size() - 1; ++j) fout << "\"" << index_paths[i][j] << "\", ";
    fout << "\"" << index_paths[i].back() << "\"]";
  }
  fout << "}}\n";
  fout.close();
}

void convert_paths_to_json(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths) {
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    cout << "{";
    cout << "\"source\":" << id_path.front() << ",";
    cout << "\"target\":" << id_path.back() << ",";
    cout << "\"length\":" << id_path.size() - 1 << ",";
    cout << "\"path\":[";
    for (int i = 0; i < (int)id_path.size() - 1; ++i) cout << id_path[i] << ",";
    cout << id_path.back() << "]";
    cout << "}" << endl;
  }
}

void convert_paths_to_json(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    fout << "{";
    fout << "\"source\":" << id_path.front() << ",";
    fout << "\"target\":" << id_path.back() << ",";
    fout << "\"length\":" << id_path.size() - 1 << ",";
    fout << "\"path\":[";
    for (int i = 0; i < (int)id_path.size() - 1; ++i) fout << id_path[i] << ",";
    fout << id_path.back() << "]";
    fout << "}" << endl;
  }
  fout.close();
}

void convert_paths_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths) {
  cout << "source,target,length,path" << endl;
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    cout << "\"" << id_path.front() << "\"";
    cout << ",\"" << id_path.back() << "\"";
    cout << ",\"" << id_path.size() - 1 << "\"";
    cout << ",\"[";
    for (int i = 0; i < (int)id_path.size() - 1; ++i) cout << id_path[i] << ", ";
    cout << id_path.back() << "]\"" << endl;
  }
}

void convert_paths_length_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths) {
  cout << "source,target,length" << endl;
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    cout << "\"" << id_path.front() << "\"";
    cout << ",\"" << id_path.back() << "\"";
    cout << ",\"" << id_path.size() - 1 << "\"";
    cout << endl;
  }
}

void convert_paths_length_col_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  fout << "length" << endl;
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    fout << "\"" << id_path.size() - 1 << "\"";
    fout << endl;
  }
  fout.close();
}

void convert_paths_length_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  fout << "source,target,length" << endl;
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    fout << "\"" << id_path.front() << "\"";
    fout << ",\"" << id_path.back() << "\"";
    fout << ",\"" << id_path.size() - 1 << "\"";
    fout << endl;
  }
  fout.close();
}

void convert_paths_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  fout << "source,target,length,path" << endl;
  for (auto &index_path : index_paths) {
    vector<long unsigned> id_path = convert_index_path_to_id_path(id_mat, index_path);
    fout << "\"" << id_path.front() << "\"";
    fout << ",\"" << id_path.back() << "\"";
    fout << ",\"" << id_path.size() - 1 << "\"";
    fout << ",\"[";
    for (int i = 0; i < (int)id_path.size() - 1; ++i) fout << id_path[i] << ", ";
    fout << id_path.back() << "]\"" << endl;
  }
  fout.close();
}

void convert_paths_to_gml(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths) {
  cout << "graph [" << endl;
  int n = id_mat.n_rows;
  for (int i = 0; i < n; ++i) {
    cout << "  node [" << endl;
    cout << "    id " << i << endl;
    cout << "    label \"" << i << "\"" << endl;
    cout << "    userid " << id_mat(i, 1) << endl;
    cout << "  ]" << endl;
  }

  unordered_set<string> edges;
  for (auto &index_path : index_paths) {
    int n = index_path.size();
    for (int i = 0; i < n - 1; ++i) {
      string edge1 = to_string(index_path[i]) + "," + to_string(index_path[i + 1]);
      string edge2 = to_string(index_path[i + 1]) + "," + to_string(index_path[i]);
      if (edges.find(edge1) != edges.end() or edges.find(edge2) != edges.end()) continue;
      cout << "  edge [" << endl;
      cout << "    source " << index_path[i] << endl;
      cout << "    target " << index_path[i + 1] << endl;
      cout << "  ]" << endl;
      edges.insert(edge1);
      edges.insert(edge2);
    }
  }
  cout << "]" << endl;
}

void convert_paths_to_gml(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename) {
  ofstream fout(filename);
  fout << "graph [" << endl;
  int n = id_mat.n_rows;
  for (int i = 0; i < n; ++i) {
    fout << "  node [" << endl;
    fout << "    id " << i << endl;
    fout << "    label \"" << i << "\"" << endl;
    fout << "    userid " << id_mat(i, 1) << endl;
    fout << "  ]" << endl;
  }

  unordered_set<string> edges;
  for (auto &index_path : index_paths) {
    int n = index_path.size();
    for (int i = 0; i < n - 1; ++i) {
      string edge1 = to_string(index_path[i]) + "," + to_string(index_path[i + 1]);
      string edge2 = to_string(index_path[i + 1]) + "," + to_string(index_path[i]);
      if (edges.find(edge1) != edges.end() or edges.find(edge2) != edges.end()) continue;
      fout << "  edge [" << endl;
      fout << "    source " << index_path[i] << endl;
      fout << "    target " << index_path[i + 1] << endl;
      fout << "  ]" << endl;
      edges.insert(edge1);
      edges.insert(edge2);
    }
  }
  fout << "]" << endl;
  fout.close();
}

void print_shortest_path(vector<int> &index_path, Mat<long unsigned> &id_mat) {
  printf("%4d: %5d (%12lu)\n", 0, index_path[0], id_mat(index_path[0], 1));
  for (int i = 1; i < (int)index_path.size(); ++i) {
    printf("%4d: %5d (%12lu)\n", i, index_path[i], id_mat(index_path[i], 1));
  }
  printf("shortest path length: %d\n", (int)index_path.size() - 1);

  printf("[%lu", id_mat(index_path[0], 1));
  for (int i = 1; i < (int)index_path.size(); ++i) {
    printf(", %lu", id_mat(index_path[i], 1));
  }
  printf("]\n");
}

void print_shortest_path(vector<int> &index_path, vector<long unsigned> &id_path) {
  printf("%4d: %5d (%12lu)\n", 0, index_path[0], id_path[0]);
  for (int i = 1; i < (int)index_path.size(); ++i) {
    printf("%4d: %5d (%12lu)\n", i, index_path[i], id_path[i]);
  }
  printf("shortest path length: %d\n", (int)index_path.size() - 1);

  printf("[%lu", id_path[0]);
  for (int i = 1; i < (int)index_path.size(); ++i) {
    printf(", %lu", id_path[i]);
  }
  printf("]\n");
}

void test_get_shortest_path_by_index(string suc_file, string id_file, int from, int to) {
  Mat<int> S;
  Mat<long unsigned> id_mat;
  printf("Reading from %s...\n", suc_file.c_str());
  S.load(suc_file, csv_ascii);
  printf("Reading from %s...\n", id_file.c_str());
  id_mat.load(id_file, csv_ascii);

  printf("Get shortest path from %d (%lu) to %d (%lu)\n", from, id_mat(from, 1), to, id_mat(to, 1));
  vector<int> index_path = get_shortest_path_by_index(S, id_mat, from, to);
  print_shortest_path(index_path, id_mat);
}

void test_get_shortest_path_by_id(string suc_file, string id_file, long unsigned id_from, long unsigned id_to) {
  Mat<int> S;
  Mat<long unsigned> id_mat;
  printf("Reading from %s...\n", suc_file.c_str());
  S.load(suc_file, csv_ascii);
  printf("Reading from %s...\n", id_file.c_str());
  id_mat.load(id_file, csv_ascii);

  pair<vector<int>, vector<long unsigned>> ret = get_shortest_path_by_id(S, id_mat, id_from, id_to);
  print_shortest_path(ret.first, ret.second);
}

void test_find_all_shortest_index_paths(string suc_file, string id_file) {
  Mat<int> S;
  Mat<long unsigned> id_mat;
  printf("Reading from %s...\n", suc_file.c_str());
  S.load(suc_file, csv_ascii);
  printf("Reading from %s...\n", id_file.c_str());
  id_mat.load(id_file, csv_ascii);
  vector<vector<int>> index_paths = find_all_shortest_index_paths(S, id_mat);

  string shortest_path_gml_file(suc_file);
  string shortest_path_csv_file(suc_file);
  string shortest_path_length_csv_file(suc_file);
  string shortest_path_length_col_csv_file(suc_file);
  string shortest_path_json_file(suc_file);
  string shortest_path_json_file_for_gml(suc_file);
  shortest_path_gml_file.replace(shortest_path_gml_file.end() - 20, shortest_path_gml_file.end(), "shortest_paths.gml");
  shortest_path_csv_file.replace(shortest_path_csv_file.end() - 20, shortest_path_csv_file.end(), "shortest_paths.csv");
  shortest_path_length_csv_file.replace(shortest_path_length_csv_file.end() - 20, shortest_path_length_csv_file.end(), "shortest_paths_length.csv");
  shortest_path_length_col_csv_file.replace(shortest_path_length_col_csv_file.end() - 20, shortest_path_length_col_csv_file.end(), "shortest_paths_length_col.csv");
  shortest_path_json_file.replace(shortest_path_json_file.end() - 20, shortest_path_json_file.end(), "shortest_paths.json");
  shortest_path_json_file_for_gml.replace(shortest_path_json_file_for_gml.end() - 20, shortest_path_json_file_for_gml.end(), "shortest_paths_for_gml.json");

  cout << endl << "Writing to " << shortest_path_gml_file << "..." << endl;
  convert_paths_to_gml(id_mat, index_paths, shortest_path_gml_file);
  cout << endl << "Writing to " << shortest_path_json_file << "..." << endl;
  convert_paths_to_json(id_mat, index_paths, shortest_path_json_file);
  cout << endl << "Writing to " << shortest_path_json_file_for_gml << "..." << endl;
  convert_paths_to_json(index_paths, shortest_path_json_file_for_gml);
  cout << endl << "Writing to " << shortest_path_csv_file << "..." << endl;
  convert_paths_to_csv(id_mat, index_paths, shortest_path_csv_file);
  cout << endl << "Writing to " << shortest_path_length_csv_file << "..." << endl;
  convert_paths_length_to_csv(id_mat, index_paths, shortest_path_length_csv_file);
  cout << endl << "Writing to " << shortest_path_length_csv_file << "..." << endl;
  convert_paths_length_col_to_csv(id_mat, index_paths, shortest_path_length_col_csv_file);
}
