(ns stats.fastica
  "An implementation of the Fast ICA algorithm (by Hyvarinen et. al.) for independent component analysis."
  (:use
     incanter.core
     incanter.stats
     [clojure.contrib.seq-utils :only (flatten)]))

;
; Convenience functions.
;

(defn mean-vector
  "Returns the mean vector of the rows of X. Assumes:
  Each row is a variable x1, x2,..., xm.
  Each column is an observation (e.g. a time-series) x1(t1), x1(t2),..., x1(tn).
  If X is an m by n matrix, a mean column vector is returned, calculated for each row, i.e. it has size m by 1."
  [X]
  (let [m (for [i (range (nrow X))]
            (mean ($ i :all X)))]
    (matrix m)))

(defn mmax
  "Finds the maximum value element of the matrix or vector x"
  [x]
  (reduce max (flatten x)))

(defn mmin
  "Finds the minimum value element of the matrix mat"
  [x]
  (reduce min (flatten x)))

(defn normalization-minmax
  "Performs min-max normalization on the value v or the matrix mat."
  [v old-min-v old-max-v new-min-v new-max-v]
  (+ 
    (* 
      (/
        (- v old-min-v)
        (- old-max-v old-min-v))
      (- new-max-v new-min-v))
    new-min-v))

(defn random-matrix
  "Creates an m by n dimensional random matrix.
  Two extra values can be supplied to specify a minimum and maximum value to normalize the matrix elements to."
  ([m n]
  (let [x (for [i (range (* m n))]
            (rand))]
    (matrix x n)))
  ([m n min-value max-value]
   (let [x (for [i (range (* m n))]
             (normalization-minmax (rand) 0 1 min-value max-value))]
     (matrix x n))))

(defn norm-euclidean
  "Calculates the Euclidean norm of a vector x."
  [x]
  (let [s (sum-of-squares x)]
    (sqrt s)))

;
; Preprocessing functions
; Centering
;

(defn center-matrix
  "Centers each row of X to achieve zero-mean for each row. Assumes:
  Each row is a variable x1, x2,..., xm.
  Each column is an observation (e.g. a time-series) x1(t1), x1(t2),..., x1(tn).
  Returns a map where :mean-vector is the mean vector of X
  and :centered-matrix is the centered matrix."
  [X]
  (let [m (mean-vector X)
        n (for [i (range (nrow X))]
            (let [row ($ i :all X)
                  mean (nth m i)]
              (minus row mean)))]
    {:mean-vector m :centered-matrix (matrix n)}))

;
; Whitening
;

(defn whiten-matrix
  "Whitens the matrix X, returning a new matrix with components that are uncorrelated. Whitening turns the covariance matrix of X into the identity matrix. Uses the eigenvalue decomposition. Assumes:
  Each row is a variable x1, x2,..., xm.
  Each column is an observation (e.g. a time-series) x1(t1), x1(t2),..., x1(tn)"
  [X]
  (let [xc (covariance (trans X)) ; X is transposed because covariance treats columns as variables.
        evd (decomp-eigenvalue xc)
        d (diag (pow (:values evd) -1/2))
        e (:vectors evd)
        wm (mmult d (trans e))
        dwm (mmult e d)
        wx (mmult wm X)]
    {:whitened-matrix wx :whitening-matrix wm :dewhitening-matrix dwm}))

;
; Contrast functions used by the negentropy calculation along with their derivatives.
; These functions are based on the recommendations in
; Independent Component Analysis: Algorithms and Applications
; by Hyvarinen and Oja.
;

(defstruct contrast-function :function :derivative)

; Power 3. Useful for sub-gaussian distribution with no outliers.
(def pow3 (struct contrast-function
                  #(* % % %)
                  #(* 3 (* % %))))

; Tanh with a = 1. Good for general purpose ICA
(def tanh (struct contrast-function
                  (fn [u] (let [a 1] (Math/tanh (* a u))))
                  (fn [u] (let [a 1
                                th (Math/tanh (* a u))]
                              (* a (- 1 (* th th)))))))

;
; Helper functions for the ICA algorithm
;

(defn weight-update-rule
  "Updates the weight (column) vector w, based on the supplied contrast function cf. Assumed that X is the whitened matrix of the mixed data where rows are mixed signals x1, x2,...,xm and columns are observations xi(t1),xi(t2),...,xi(tn).
  Inspired by the Java implementation of the FastICA algorithm which uses the stabilized version of FastICA."
  [X w cf]
  (let [wx (mmult (trans w) X) ; wx, row vector
        gwx (map (:function cf) wx)   ; g(wx), column vector
        gdwx (map (:derivative cf) wx) ; g'(wx), column vector
        beta (div (mmult wx gwx) (ncol X)) ; single value
        exgwx (div (mmult X gwx) (ncol X)) ; E{x g(w'x)}, column vector
        egdwx (div (sum gdwx) (ncol X))    ; E{g'(w'x)}, single value
        new-w (minus
                w
                (div
                  (minus exgwx (mult beta w))
                  (- egdwx beta)))]
    new-w))

(defn fastica-one-unit
  "Finds a weight vector that will transform X into one independent component.
  Uses the contrast function cf.
  Uses deflation decorrelation according to the matrix B. Assumes that X is centered and whitened. The weight vector is returned when the maximum number of iterations are reached (max-it) or the weight vector has converged with an error less than eps."
  [X B cf max-it eps]
  (loop [cur-it 0
         w (random-matrix (nrow X) 1)
         w-old (matrix 0 (nrow X) 1)]
    ; Decorrelate the weight vector with the previous vectors
    ; Deflation method
    (let [w (minus w (mmult (mmult B (trans B)) w))
          w (div w (norm-euclidean w))]
      (println "Iteration" cur-it)
      ; Check for convergence or maximum iterations
      (if (or (>= cur-it max-it)
              (or
                (< (norm-euclidean (minus w w-old)) eps)
                (< (norm-euclidean (plus w w-old)) eps)))
        w
        (let [new-w (weight-update-rule X w cf)]
          (recur (inc cur-it) (div new-w (norm-euclidean new-w)) w))))))

(defn fastica
  "The FastICA algoritm. Takes a data matrix, X, consisting of mixed signals, x1, x2,..., xm, and tries to approximate the original signals s1, s2,..., sm.
  The algorithm uses an approximation to negentropy with the given derivative and double derivative (cf) of the contrast function G.
  It is assumed that the rows of X are variables and the columns are the observations.
  Currently assumes number of rows = number of components to find.
  Currently uses the deflation approach for de-correlating the found weight vectors."
  [X & options]
  (let [; Begin options
        opts (apply hash-map options)
        eps (if (:eps opts) (:eps opts) 0.0001)
        max-it (if (:max-it opts) (:max-it opts) 100)
        cf (if (:cf opts) (:cf opts) pow3)
        ; End options
        ; Begin pre-processing
        num-ica (nrow X)
        centering (center-matrix X)
        whitening (whiten-matrix (:centered-matrix centering))
        X-whitened (:whitened-matrix whitening)]
        ; End pre-processing
    ; Main loop, finds independent components one at a time
    (loop [cur-ic 0
           B (matrix 0 num-ica num-ica)  ; B = Decorrelation matrix
           W (matrix 0 num-ica num-ica)] ; W = final weight matrix
      (println "Finding component" cur-ic)
      (if (= cur-ic num-ica) ; All independent components are found
        W ; Return the final de-mixing matrix
        (let [w (fastica-one-unit X-whitened B cf max-it eps)]
          (recur
            (inc cur-ic)
            (bind-columns ($ :all (range 0 cur-ic) B)
                          w
                          ($ :all (range (inc cur-ic) num-ica) B))
            (bind-rows ($ (range 0 cur-ic) :all W)
                       (mmult (trans w) (:whitening-matrix whitening))
                       ($ (range (inc cur-ic) num-ica) :all W))))))))
