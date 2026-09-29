module float_type = f64
type t = float_type.t

let pi: t = 3.14159265358979323846

-- Calculation method
let meth_gauss_seidel: i64 = 1
let meth_jacobi: i64 = 2

-- Perturbation function
let func_f0: i64 = 1
let func_fpisin: i64 = 2

-- Termination condition
let term_acc: i64 = 1
let term_iter: i64 = 2

let matrix_size (interlines: i64): i64 = interlines * 8 + 8

let get_residuum [m] (old: [m][m]t) (new: [m][m]t): t =
  map2 (\row_old row_new -> map2 (\a b -> float_type.abs (a - b)) row_old row_new) old new
  |> flatten
  |> float_type.maximum

let jacobi_sweep [m] (h: t) (func: i64) (matrix: [m][m]t): (*[m][m]t, t) =
  let n = m - 1
  let pih = pi * h
  let fpisin = 0.25 * 2.0 * pi * pi * h * h
  let new =
    tabulate_2d m m
      (\i j ->
         if i == 0 || i == n || j == 0 || j == n
         then matrix[i][j]
         else let star = 0.25 * (matrix[i - 1][j] + matrix[i][j - 1] +
                                  matrix[i][j + 1] + matrix[i + 1][j])
              in if func == func_fpisin
                 then star + fpisin * float_type.sin (pih * float_type.i64 i) * float_type.sin (pih * float_type.i64 j)
                 else star)
  in (new, get_residuum matrix new)

let gauss_seidel_sweep [m] (h: t) (func: i64) (matrix0: *[m][m]t): (*[m][m]t, t) =
  let n = m - 1
  let pih = pi * h
  let fpisin = 0.25 * 2.0 * pi * pi * h * h
  in loop (matrix, maxres) = (matrix0, 0.0) for i_ < (n - 1) do
       let i = i_ + 1
       in loop (matrix, maxres) = (matrix, maxres) for j_ < (n - 1) do
            let j = j_ + 1
            let star = 0.25 * (matrix[i - 1, j] + matrix[i, j - 1] +
                                matrix[i, j + 1] + matrix[i + 1, j])
            let star = if func == func_fpisin
                       then star + fpisin * float_type.sin (pih * float_type.i64 i) * float_type.sin (pih * float_type.i64 j)
                       else star
            let residuum = float_type.abs (matrix[i, j] - star)
            let matrix[i, j] = star
            in (matrix, float_type.max maxres residuum)

let sweep [m] (h: t) (func: i64) (method: i64) (matrix: *[m][m]t): (*[m][m]t, t) =
  if method == meth_jacobi
  then jacobi_sweep h func matrix
  else gauss_seidel_sweep h func matrix

entry init_matrices (interlines: i64) (func: i64): [][]t =
  let n = matrix_size interlines
  let m = n + 1
  let h = 1.0 / float_type.i64 n
  in tabulate_2d m m
       (\i j ->
          if func != func_f0
          then 0.0
          else if i == 0 && j == 0
          then 1.0
          else if (i == n && j == 0) || (i == 0 && j == n)
          then 0.0
          else if j == 0
          then 1.0 - h * float_type.i64 i
          else if j == n
          then h * float_type.i64 i
          else if i == 0
          then 1.0 - h * float_type.i64 j
          else if i == n
          then h * float_type.i64 j
          else 0.0)

entry calculate (method: i64) (func: i64) (term: i64) (acc_iter: t) (matrix0: *[][]t): ([][]t, t, i64) =
  assert (method == meth_gauss_seidel || method == meth_jacobi) (
  assert (func == func_f0 || func == func_fpisin) (
  assert (term == term_acc || term == term_iter) (
  let m = length matrix0
  let n = m - 1
  let h = 1.0 / float_type.i64 n
  let term_iteration_max = float_type.to_i64 acc_iter
  in loop (matrix, res, iter) = (matrix0, float_type.inf, 0)
     while (iter == 0 ||
            (if term == term_iter then iter < term_iteration_max else res >= acc_iter)) do
       let (matrix', res') = sweep h func method matrix
       in (matrix', res', iter + 1))))
