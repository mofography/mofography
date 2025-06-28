use ndarray::Array3;
use numpy::{PyArray3, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::prelude::*;
use hashbrown::HashMap;
use numpy::{PyArray1, PyArray2, PyReadonlyArray3};
use pyo3::types::{PyList, PyTuple};
use ordered_float::NotNan;
use parking_lot::Mutex;
use std::sync::Arc;

/// Compute a distance grid (minimum atom surface distance at each grid point) for a periodic
/// structure.
///
/// Parameters
/// ----------
/// fpoints : numpy.ndarray[float64] shape=(n_atoms, 3)
///     Fractional coordinates of atoms shifted so that the unit cell spans [-0.5, 0.5] in each
///     dimension (i.e. the same convention used in the original Python implementation).
/// cell : numpy.ndarray[float64] shape=(3, 3)
///     3×3 matrix whose columns are the crystallographic lattice vectors **in Å**.
/// radii : numpy.ndarray[float64] shape=(n_atoms,)
///     van-der-Waals radii (or any other probe radius) for each atom **in Å**.
/// nx, ny, nz : int
///     Number of grid samples in the a, b, and c directions respectively (identical to the
///     `np.linspace(-0.5, 0.5, n?)` convention in the Python reference).
///
/// Returns
/// -------
/// numpy.ndarray[float32] shape=(nx, ny, nz)
///     Distance grid with the same semantics as the reference Python version (float32 to save
///     memory).
#[pyfunction]
#[pyo3(text_signature = "(fpoints, cell, radii, nx, ny, nz)")]
fn dgrid_from_atoms(
    py: Python<'_>,
    fpoints: PyReadonlyArray2<'_, f64>,
    cell: PyReadonlyArray2<'_, f64>,
    radii: PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
) -> PyResult<Py<PyArray3<f32>>> {
    // Convert input arrays to Rust-native views.
    let fpoints = fpoints.as_array(); // Shape (n_atoms, 3)
    let cell = cell.as_array();       // Shape (3, 3)
    let radii = radii.as_array();     // Shape (n_atoms,)

    let total_points = nx * ny * nz;

    // Pre-compute cell components for fast access.
    let c00 = cell[[0, 0]];
    let c01 = cell[[0, 1]];
    let c02 = cell[[0, 2]];
    let c10 = cell[[1, 0]];
    let c11 = cell[[1, 1]];
    let c12 = cell[[1, 2]];
    let c20 = cell[[2, 0]];
    let c21 = cell[[2, 1]];
    let c22 = cell[[2, 2]];

    // Parallel computation across all grid points.
    let buffer: Vec<f32> = (0..total_points)
        .into_par_iter()
        .map(|idx| {
            // Unravel flat index into 3-D indices (Fortran order).
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            // Fractional coordinates on the uniform grid (inclusive end-points like linspace).
            let x_frac = -0.5 + (i as f64) / ((nx - 1) as f64);
            let y_frac = -0.5 + (j as f64) / ((ny - 1) as f64);
            let z_frac = -0.5 + (k as f64) / ((nz - 1) as f64);

            // Track the minimum atom-surface distance.
            let mut min_dist = f64::MAX;

            for (atom_idx, atom_frac) in fpoints.outer_iter().enumerate() {
                // Apply minimum-image convention in fractional space.
                let mut dx = x_frac - atom_frac[0];
                let mut dy = y_frac - atom_frac[1];
                let mut dz = z_frac - atom_frac[2];

                dx -= dx.round();
                dy -= dy.round();
                dz -= dz.round();

                // Convert to Cartesian coordinates: r = cell · diff
                let cart_x = c00 * dx + c01 * dy + c02 * dz;
                let cart_y = c10 * dx + c11 * dy + c12 * dz;
                let cart_z = c20 * dx + c21 * dy + c22 * dz;

                let dist = (cart_x * cart_x + cart_y * cart_y + cart_z * cart_z).sqrt()
                    - radii[atom_idx];

                if dist < min_dist {
                    min_dist = dist;
                }
            }

            // Cast to f32 – this is enough precision for the grid and halves memory usage.
            min_dist as f32
        })
        .collect();

    // Reshape the flat Vec into a 3-D ndarray (row-major as in NumPy).
    let array = Array3::from_shape_vec((nx, ny, nz), buffer)
        .expect("Shape mismatch while reshaping the output distance grid");

    // Convert to a NumPy array owned by Python and return.
    Ok(PyArray3::from_owned_array(py, array).to_owned())
}

#[pyfunction]
#[pyo3(text_signature = "(fpoints, cell, epsilons, sigmas, nx, ny, nz, cutoff)")]
fn egrid_from_atoms(
    py: Python<'_>,
    fpoints: PyReadonlyArray2<'_, f64>,
    cell: PyReadonlyArray2<'_, f64>,
    epsilons: PyReadonlyArray1<'_, f64>,
    sigmas: PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    cutoff: f64,
) -> PyResult<Py<PyArray3<f64>>> {
    let fpoints = fpoints.as_array();
    let cell = cell.as_array();
    let epsilons = epsilons.as_array();
    let sigmas = sigmas.as_array();

    let total_points = nx * ny * nz;
    let cutoff_sq = cutoff * cutoff;

    // Pre-compute cell components for fast access.
    let c00 = cell[[0, 0]];
    let c01 = cell[[0, 1]];
    let c02 = cell[[0, 2]];
    let c10 = cell[[1, 0]];
    let c11 = cell[[1, 1]];
    let c12 = cell[[1, 2]];
    let c20 = cell[[2, 0]];
    let c21 = cell[[2, 1]];
    let c22 = cell[[2, 2]];

    // Parallel computation across all grid points.
    let buffer: Vec<f64> = (0..total_points)
        .into_par_iter()
        .map(|idx| {
            // Unravel flat index into 3-D indices (Fortran order).
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            // Fractional coordinates on the uniform grid (inclusive end-points like linspace).
            let x_frac = -0.5 + (i as f64) / ((nx - 1) as f64);
            let y_frac = -0.5 + (j as f64) / ((ny - 1) as f64);
            let z_frac = -0.5 + (k as f64) / ((nz - 1) as f64);

            let mut energy = 0.0_f64;

            for (atom_idx, atom_frac) in fpoints.outer_iter().enumerate() {
                // Minimum-image convention in fractional space.
                let mut dx = x_frac - atom_frac[0];
                let mut dy = y_frac - atom_frac[1];
                let mut dz = z_frac - atom_frac[2];

                dx -= dx.round();
                dy -= dy.round();
                dz -= dz.round();

                // Convert to Cartesian coordinates.
                let cart_x = c00 * dx + c01 * dy + c02 * dz;
                let cart_y = c10 * dx + c11 * dy + c12 * dz;
                let cart_z = c20 * dx + c21 * dy + c22 * dz;

                let r_sq = cart_x * cart_x + cart_y * cart_y + cart_z * cart_z;

                // Skip if beyond cutoff.
                if r_sq > cutoff_sq {
                    continue;
                }

                // Avoid singularity when r -> 0.
                let r = if r_sq < 1e-24 { 1e-12 } else { r_sq.sqrt() };

                let sigma = sigmas[atom_idx];
                let epsilon = epsilons[atom_idx];

                let sr = sigma / r;
                let sr6 = sr.powi(6);
                // Lennard-Jones 12-6 potential.
                energy += 4.0 * epsilon * (sr6 * sr6 - sr6);
            }

            energy
        })
        .collect();

    let array = Array3::from_shape_vec((nx, ny, nz), buffer)
        .expect("Shape mismatch while reshaping the output energy grid");

    Ok(PyArray3::from_owned_array(py, array).to_owned())
}

#[derive(Debug, Clone)]
struct WindowPoint {
    coord: [f64; 3],
    region: u32,
    distance: f32,
}

// Simple union-find for clustering
#[derive(Debug)]
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        
        if root_x != root_y {
            if self.rank[root_x] < self.rank[root_y] {
                self.parent[root_x] = root_y;
            } else if self.rank[root_x] > self.rank[root_y] {
                self.parent[root_y] = root_x;
            } else {
                self.parent[root_y] = root_x;
                self.rank[root_x] += 1;
            }
        }
    }

    fn get_clusters(&mut self) -> HashMap<usize, Vec<usize>> {
        let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..self.parent.len() {
            let root = self.find(i);
            clusters.entry(root).or_default().push(i);
        }
        clusters
    }
}

fn euclidean_distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn cluster_points(points: &[WindowPoint], distance_threshold: f64) -> Vec<Vec<usize>> {
    if points.is_empty() {
        return vec![];
    }
    
    let n = points.len();
    let mut uf = UnionFind::new(n);
    
    // Connect points within distance threshold
    for i in 0..n {
        for j in (i + 1)..n {
            if euclidean_distance(&points[i].coord, &points[j].coord) < distance_threshold {
                uf.union(i, j);
            }
        }
    }
    
    // Extract clusters
    let cluster_map = uf.get_clusters();
    cluster_map.into_values().collect()
}

/// Apply a periodic (wrap-around) Gaussian filter to a 3D array.
///
/// This is a memory-efficient implementation that avoids creating multiple copies
/// of the input array, unlike the Python version that uses np.roll.
///
/// Parameters
/// ----------
/// input : numpy.ndarray[float64] shape=(nx, ny, nz)
///     Input 3D array to be filtered
/// sigma : float
///     Standard deviation of the Gaussian kernel
///
/// Returns
/// -------
/// numpy.ndarray[float64] shape=(nx, ny, nz)
///     Filtered array with same shape as input
#[pyfunction]
#[pyo3(text_signature = "(input, sigma)")]
fn periodic_gaussian_filter(
    py: Python<'_>,
    input: PyReadonlyArray3<'_, f64>,
    sigma: f64,
) -> PyResult<Py<PyArray3<f64>>> {
    use ndarray::{Axis, s};

    let inp = input.as_array();
    let (nx, ny, nz) = inp.dim();

    // Build 1-D kernel (truncate at 3σ)
    let radius = (3.0 * sigma).ceil() as i32;
    let mut kernel: Vec<f64> = (-radius..=radius)
        .map(|k| {
            let kf = k as f64;
            (-0.5 * kf * kf / (sigma * sigma)).exp()
        })
        .collect();
    let norm: f64 = kernel.iter().sum();
    for w in &mut kernel {
        *w /= norm;
    }

    // Temp buffers
    let mut tmp_x = Array3::<f64>::zeros((nx, ny, nz));
    let mut tmp_y = Array3::<f64>::zeros((nx, ny, nz));

    // ---- Pass 1 : X axis ----
    tmp_x
        .axis_iter_mut(Axis(1)) // iterate over y
        .enumerate()
        .for_each(|(y_idx, mut slab)| {
            // slab has shape (nx, nz) for fixed y
            for z in 0..nz {
                for x in 0..nx {
                    let mut acc = 0.0;
                    for (offset, &w) in (-radius..=radius).zip(kernel.iter()) {
                        let xi = (x as i32 + offset).rem_euclid(nx as i32) as usize;
                        acc += w * inp[[xi, y_idx, z]];
                    }
                    slab[[x, z]] = acc;
                }
            }
        });

    // ---- Pass 2 : Y axis ----
    tmp_y
        .axis_iter_mut(Axis(0)) // iterate over x
        .enumerate()
        .for_each(|(x_idx, mut slab)| {
            // slab shape (ny, nz) for fixed x
            for z in 0..nz {
                for y in 0..ny {
                    let mut acc = 0.0;
                    for (offset, &w) in (-radius..=radius).zip(kernel.iter()) {
                        let yj = (y as i32 + offset).rem_euclid(ny as i32) as usize;
                        acc += w * tmp_x[[x_idx, yj, z]];
                    }
                    slab[[y, z]] = acc;
                }
            }
        });

    // ---- Pass 3 : Z axis ----
    let mut out = Array3::<f64>::zeros((nx, ny, nz));
    out
        .axis_iter_mut(Axis(2)) // iterate over z
        .enumerate()
        .for_each(|(z_idx, mut slab)| {
            // slab shape (nx, ny) for fixed z
            for y in 0..ny {
                for x in 0..nx {
                    let mut acc = 0.0;
                    for (offset, &w) in (-radius..=radius).zip(kernel.iter()) {
                        let zk = (z_idx as i32 + offset).rem_euclid(nz as i32) as usize;
                        acc += w * tmp_y[[x, y, zk]];
                    }
                    slab[[x, y]] = acc;
                }
            }
        });

    Ok(PyArray3::from_owned_array(py, out).to_owned())
}

#[pyfunction]
#[pyo3(text_signature = "(labels, dgrid, padding_thickness, cell, distance_threshold)")]
fn window_connections(
    py: Python<'_>,
    labels: PyReadonlyArray3<'_, u32>,
    dgrid: PyReadonlyArray3<'_, f32>,
    padding_thickness: usize,
    cell: PyReadonlyArray2<'_, f64>,
    distance_threshold: f64,
) -> PyResult<PyObject> {
    // Convert to ndarray views
    let labels = labels.as_array();
    let dgrid = dgrid.as_array();
    let cell = cell.as_array();

    let (nx, ny, nz) = (labels.dim().0, labels.dim().1, labels.dim().2);

    // Helper closure to transform voxel indices -> Cartesian coordinate
    let to_cart = |i: usize, j: usize, k: usize| -> [f64; 3] {
        // Convert padded grid indices to fractional coordinates relative to original grid
        let frac_x = (i as f64 - padding_thickness as f64) / (nx as f64 - 2.0 * padding_thickness as f64);
        let frac_y = (j as f64 - padding_thickness as f64) / (ny as f64 - 2.0 * padding_thickness as f64);
        let frac_z = (k as f64 - padding_thickness as f64) / (nz as f64 - 2.0 * padding_thickness as f64);
        
        // Convert to Cartesian coordinates using cell matrix
        let cx = cell[[0, 0]] * frac_x + cell[[0, 1]] * frac_y + cell[[0, 2]] * frac_z;
        let cy = cell[[1, 0]] * frac_x + cell[[1, 1]] * frac_y + cell[[1, 2]] * frac_z;
        let cz = cell[[2, 0]] * frac_x + cell[[2, 1]] * frac_y + cell[[2, 2]] * frac_z;
        [cx, cy, cz]
    };

    // Map edge key -> Vec of boundary points
    let mut edge_map: HashMap<(u32, u32), Vec<WindowPoint>> = HashMap::new();

    // Scan grid for boundary voxels between different regions
    // Check all 6 neighbors to find region boundaries
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let label_ijk = labels[[i, j, k]];
                if label_ijk == 0 {
                    continue;
                }
                
                let neighbors = [
                    (i.wrapping_sub(1), j, k), // -x
                    (i + 1, j, k),              // +x
                    (i, j.wrapping_sub(1), k), // -y
                    (i, j + 1, k),              // +y
                    (i, j, k.wrapping_sub(1)), // -z
                    (i, j, k + 1),              // +z
                ];
                
                for (ni, nj, nk) in neighbors {
                    if ni < nx && nj < ny && nk < nz {
                        let neighbor_label = labels[[ni, nj, nk]];
                        if neighbor_label != 0 && neighbor_label != label_ijk {
                            // Found boundary between two regions
                            let (reg1, reg2) = if label_ijk < neighbor_label {
                                (label_ijk, neighbor_label)
                            } else {
                                (neighbor_label, label_ijk)
                            };
                            
                            let coord = to_cart(i, j, k);
                            let point = WindowPoint {
                                coord,
                                region: label_ijk,
                                distance: dgrid[[i, j, k]],
                            };
                            
                            edge_map.entry((reg1, reg2)).or_default().push(point);
                        }
                    }
                }
            }
        }
    }

    // Process each edge in parallel and build output
    let connections: Vec<(u32, u32, Vec<[f64; 3]>, Vec<f64>)> = edge_map
        .into_iter()
        .collect::<Vec<_>>()
        .into_par_iter()
        .filter_map(|((reg1, reg2), boundary_points)| {
            if boundary_points.is_empty() {
                return None;
            }

            // Cluster boundary points using distance threshold
            let clusters = cluster_points(&boundary_points, distance_threshold);
            
            if clusters.is_empty() {
                return None;
            }

            let mut window_centers = Vec::new();
            let mut window_radii = Vec::new();

            for cluster_indices in clusters {
                if cluster_indices.len() < 2 {
                    continue; // Need points from both regions
                }

                // Separate points by region within this cluster
                let mut reg1_points = Vec::new();
                let mut reg2_points = Vec::new();
                
                for &idx in &cluster_indices {
                    let point = &boundary_points[idx];
                    if point.region == reg1 {
                        reg1_points.push(point);
                    } else if point.region == reg2 {
                        reg2_points.push(point);
                    }
                }

                if reg1_points.is_empty() || reg2_points.is_empty() {
                    continue; // Need points from both regions
                }

                // Find point with maximum distance value in each region
                let max_reg1 = reg1_points
                    .iter()
                    .max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())
                    .unwrap();
                let max_reg2 = reg2_points
                    .iter()
                    .max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())
                    .unwrap();

                // Calculate window center as midpoint of the two maxima
                let center = [
                    0.5 * (max_reg1.coord[0] + max_reg2.coord[0]),
                    0.5 * (max_reg1.coord[1] + max_reg2.coord[1]),
                    0.5 * (max_reg1.coord[2] + max_reg2.coord[2]),
                ];

                // Use average of the two maximum distance values as radius
                let radius = 0.5 * (max_reg1.distance + max_reg2.distance) as f64;

                window_centers.push(center);
                window_radii.push(radius);
            }

            if window_centers.is_empty() {
                None
            } else {
                Some((reg1, reg2, window_centers, window_radii))
            }
        })
        .collect();

    // Convert to Python objects
    let py_out = PyList::empty(py);

    for (reg1, reg2, centers, radii) in connections {
        // Build numpy arrays
        let regs_arr = PyArray1::<u32>::from_vec(py, vec![reg1, reg2]);
        
        // Convert [f64; 3] to Vec<f64> for each center
        let centers_vec: Vec<Vec<f64>> = centers
            .iter()
            .map(|center| center.to_vec())
            .collect();
        let centers_arr = PyArray2::<f64>::from_vec2(py, &centers_vec)?;
        let radii_arr = PyArray1::<f64>::from_vec(py, radii);

        // Create tuple matching Python format: [regs, connected, type, centers, radii]
        let tup = PyTuple::new(
            py,
            &[
                regs_arr.into_py(py),
                true.into_py(py),
                "internal".into_py(py),
                centers_arr.into_py(py),
                radii_arr.into_py(py),
            ],
        );
        py_out.append(tup)?;
    }

    Ok(py_out.into())
}

const BUCKETS: usize = 65536;

#[pyfunction]
#[pyo3(text_signature = "(image, seeds, mask=None, n_threads=None)")]
fn seeded_watershed_bucket(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    mask: Option<PyReadonlyArray3<'_, u8>>, // optional
    n_threads: Option<usize>,
) -> PyResult<Py<PyArray3<u32>>> {
    let img = image.as_array();
    let seeds_arr = seeds.as_array();

    if img.dim() != seeds_arr.dim() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "`image` and `seeds` must have identical shapes",
        ));
    }
    if let Some(ref m) = mask {
        if m.as_array().dim() != img.dim() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "`mask` shape must match `image`",
            ));
        }
    }

    if let Some(n) = n_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    // Determine min and max of the image to scale into 0..65535
    let (mut min_val, mut max_val) = (f32::MAX, f32::MIN);
    img.iter().for_each(|&v| {
        if v < min_val {
            min_val = v;
        }
        if v > max_val {
            max_val = v;
        }
    });
    // avoid div-by-zero
    let range = if (max_val - min_val).abs() < 1e-12 {
        1.0
    } else {
        max_val - min_val
    };

    let to_bucket = |v: f32| -> usize {
        let norm = ((v - min_val) / range * (BUCKETS as f32 - 1.0)).round() as usize;
        norm
    };

    // Precompute bucket index for every voxel to avoid recalculating
    let bucket_idx_vec: Vec<usize> = img
        .iter()
        .map(|&v| to_bucket(v))
        .collect();

    let (nx, ny, nz) = img.dim();
    let total_voxels = nx * ny * nz;

    use std::sync::atomic::{AtomicU32, Ordering};

    let labels: Vec<AtomicU32> = (0..total_voxels).map(|_| AtomicU32::new(0)).collect();

    let buckets: Arc<Vec<Mutex<Vec<usize>>>> = Arc::new(
        (0..BUCKETS).map(|_| Mutex::new(Vec::<usize>::new())).collect(),
    );

    let to_flat = |i: usize, j: usize, k: usize| -> usize { i * ny * nz + j * nz + k };

    // Materialise mask (if any) into a flat Vec<bool> that is Sync/Send.
    let mask_vec: Option<Vec<bool>> = mask.as_ref().map(|m| {
        m.as_array()
            .iter()
            .map(|&v| v != 0)
            .collect::<Vec<bool>>()
    });

    // Initialise with seeds
    for ((i, j, k), &label) in seeds_arr.indexed_iter() {
        if label == 0 {
            continue;
        }
        if let Some(ref mask_v) = mask_vec {
            if !mask_v[to_flat(i, j, k)] {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Seed label found in location where mask==0",
                ));
            }
        }
        let idx = to_flat(i, j, k);
        labels[idx].store(label, Ordering::Relaxed);
        let b = bucket_idx_vec[idx];
        buckets[b].lock().push(idx);
    }

    let neighbor_offsets: &[(isize, isize, isize)] = &[
        (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
        (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1),
        (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1),                 (0, 0, 1),
        (0, 1, -1), (0, 1, 0), (0, 1, 1),
        (1, -1, -1), (1, -1, 0), (1, -1, 1),
        (1, 0, -1), (1, 0, 0), (1, 0, 1),
        (1, 1, -1), (1, 1, 0), (1, 1, 1),
    ];

    // Main flooding loop from highest bucket downwards
    for level in (0..BUCKETS).rev() {
        let frontier: Vec<usize> = {
            let mut guard = buckets[level].lock();
            std::mem::take(&mut *guard)
        };
        if frontier.is_empty() {
            continue;
        }

        frontier.par_iter().for_each(|&idx| {
            let current_label = labels[idx].load(Ordering::Relaxed);
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            for &(dx, dy, dz) in neighbor_offsets {
                let ni = i as isize + dx;
                let nj = j as isize + dy;
                let nk = k as isize + dz;
                if ni < 0 || nj < 0 || nk < 0 {
                    continue;
                }
                let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
                if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                    continue;
                }
                if let Some(ref mask_v) = mask_vec {
                    if !mask_v[to_flat(ni_u, nj_u, nk_u)] {
                        continue;
                    }
                }
                let nidx = to_flat(ni_u, nj_u, nk_u);
                if labels[nidx]
                    .compare_exchange(0, current_label, Ordering::Relaxed, Ordering::Relaxed)
                    .is_ok()
                {
                    let b = bucket_idx_vec[nidx];
                    buckets[b].lock().push(nidx);
                }
            }
        });
    }

    // Convert labels into ndarray
    let mut out_vec: Vec<u32> = vec![0u32; total_voxels];
    out_vec
        .iter_mut()
        .enumerate()
        .for_each(|(idx, v)| *v = labels[idx].load(Ordering::Relaxed));

    let out_arr = Array3::from_shape_vec((nx, ny, nz), out_vec)
        .expect("shape mismatch in output");

    Ok(PyArray3::from_owned_array(py, out_arr).to_owned())
}

// ============================================================================
// Heap-based single-threaded watershed (original version, retained for backwards
// compatibility and low-memory usage)
// ============================================================================

#[pyfunction]
#[pyo3(text_signature = "(image, seeds)")]
fn seeded_watershed(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
) -> PyResult<Py<PyArray3<u32>>> {
    heap_watershed_impl(py, image, seeds, None)
}

#[pyfunction]
#[pyo3(text_signature = "(image, seeds, mask)")]
fn seeded_watershed_masked(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    mask: PyReadonlyArray3<'_, u8>,
) -> PyResult<Py<PyArray3<u32>>> {
    heap_watershed_impl(py, image, seeds, Some(mask))
}

// Shared helper implementing the original heap algorithm ------------------------------------
fn heap_watershed_impl(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    mask: Option<PyReadonlyArray3<'_, u8>>,
) -> PyResult<Py<PyArray3<u32>>> {
    let img = image.as_array();
    let seeds_arr = seeds.as_array();
    if img.dim() != seeds_arr.dim() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "`image` and `seeds` must have identical shapes",
        ));
    }

    if let Some(ref m) = mask {
        if m.as_array().dim() != img.dim() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "`mask` shape must match `image`",
            ));
        }
    }

    let (nx, ny, nz) = img.dim();
    let total_voxels = nx * ny * nz;

    // Flattened labels
    let mut labels: Vec<u32> = vec![0; total_voxels];
    let to_flat = |i: usize, j: usize, k: usize| -> usize { i * ny * nz + j * nz + k };

    use std::collections::BinaryHeap;
    let mut heap: BinaryHeap<(NotNan<f32>, usize)> = BinaryHeap::new();

    // Seed initialisation
    for ((i, j, k), &label) in seeds_arr.indexed_iter() {
        if label == 0 {
            continue;
        }
        if let Some(ref m) = mask {
            if m.as_array()[[i, j, k]] == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Seed label found in location where mask==0",
                ));
            }
        }
        let idx = to_flat(i, j, k);
        labels[idx] = label;
        let priority = NotNan::new(img[[i, j, k]]).unwrap();
        heap.push((priority, idx));
    }

    // 26-NN offsets
    let neighbor_offsets: &[(isize, isize, isize)] = &[
        (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
        (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1),
        (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1),                 (0, 0, 1),
        (0, 1, -1), (0, 1, 0), (0, 1, 1),
        (1, -1, -1), (1, -1, 0), (1, -1, 1),
        (1, 0, -1), (1, 0, 0), (1, 0, 1),
        (1, 1, -1), (1, 1, 0), (1, 1, 1),
    ];

    while let Some((_prio, idx)) = heap.pop() {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let current_label = labels[idx];

        for &(dx, dy, dz) in neighbor_offsets {
            let ni = i as isize + dx;
            let nj = j as isize + dy;
            let nk = k as isize + dz;
            if ni < 0 || nj < 0 || nk < 0 {
                continue;
            }
            let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
            if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                continue;
            }
            if let Some(ref m) = mask {
                if m.as_array()[[ni_u, nj_u, nk_u]] == 0 {
                    continue;
                }
            }
            let nidx = to_flat(ni_u, nj_u, nk_u);
            if labels[nidx] == 0 {
                labels[nidx] = current_label;
                let priority = NotNan::new(img[[ni_u, nj_u, nk_u]]).unwrap();
                heap.push((priority, nidx));
            }
        }
    }

    let out_arr = Array3::from_shape_vec((nx, ny, nz), labels).expect("shape mismatch");
    Ok(PyArray3::from_owned_array(py, out_arr).to_owned())
}

#[pymodule]
fn mgr_rust(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // Enable better error reporting when the module is imported from Python.
    pyo3_log::init();

    m.add_function(wrap_pyfunction!(dgrid_from_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(egrid_from_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(window_connections, m)?)?;
    m.add_function(wrap_pyfunction!(periodic_gaussian_filter, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_watershed_bucket, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_watershed, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_watershed_masked, m)?)?;
    m.add("__all__", ["dgrid_from_atoms", "egrid_from_atoms", "window_connections", "periodic_gaussian_filter", "seeded_watershed_bucket", "seeded_watershed", "seeded_watershed_masked"])?;

    // By returning Ok(()), the module initialization is successful.
    Ok(())
} 