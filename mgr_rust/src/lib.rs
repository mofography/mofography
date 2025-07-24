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

#[pyfunction]
#[pyo3(signature = (grid, spacing, h, threshold_abs=None, min_distance=0.0, cell=None, refine_peaks=false, max_steps=100))]
fn detect_maxima_h_morph(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    spacing: f32,
    h: f32,
    threshold_abs: Option<f32>,
    min_distance: f32,                       // Å
    cell: Option<PyReadonlyArray2<'_, f64>>, // 3×3 matrix with unit vectors as columns
    refine_peaks: bool,                      // Whether to refine peak locations using gradient ascent
    max_steps: usize,                        // Max steps for gradient ascent
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array2, Array3};
    use rayon::prelude::*;

    // Start with original grid
    let mut img: Array3<f32> = grid.as_array().to_owned();
    let (nx, ny, nz) = img.dim();

    // Apply absolute threshold if provided
    if let Some(tabs) = threshold_abs {
        img.mapv_inplace(|v| if v < tabs { -1e6 } else { v });
    }

    let total_voxels = nx * ny * nz;

    // Early return for tiny grids
    if total_voxels < 27 {
        let empty = Array2::<usize>::zeros((0, 3));
        return Ok((
            PyArray2::from_owned_array(py, empty.clone()).to_owned(),
            PyArray2::from_owned_array(py, empty.clone()).to_owned(),
            PyArray2::from_owned_array(py, empty).to_owned(),
        ));
    }

    // ----------- Robust morphological reconstruction by dilation (iterative) ---------------------

    // Create marker image (f - h) - this is the seed for reconstruction
    let marker = img.mapv(|v| v - h);
    let mut recon = marker.clone();

    // 26-connectivity neighbors for reconstruction (no periodic wrapping)
    let neighbors: &[(isize, isize, isize)] = &[
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

    // Iterative geodesic reconstruction until convergence
    // This is more robust than queue-based for numerical precision
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        
        // Forward pass
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let current_mask = img[[i, j, k]];
                    
                    // Skip thresholded-out regions
                    if current_mask <= -1e5 {
                        continue;
                    }
                    
                    let mut max_neighbor = recon[[i, j, k]];
                    
                    // Check all 26 neighbors
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;

                        // Check bounds (no wrapping)
                        if ni < 0 || nj < 0 || nk < 0 {
                            continue;
                        }
                        let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
                        if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                            continue;
                        }
                        
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    
                    // Geodesic reconstruction: constrain by mask
                    let new_val = max_neighbor.min(current_mask);
                    
                    // Use a small epsilon for floating point comparison
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
        
        // Backward pass for better convergence
        for i in (0..nx).rev() {
            for j in (0..ny).rev() {
                for k in (0..nz).rev() {
                    let current_mask = img[[i, j, k]];
                    
                    // Skip thresholded-out regions
                    if current_mask <= -1e5 {
                        continue;
                    }
                    
                    let mut max_neighbor = recon[[i, j, k]];
                    
                    // Check all 26 neighbors
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;

                        // Check bounds (no wrapping)
                        if ni < 0 || nj < 0 || nk < 0 {
                            continue;
                        }
                        let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
                        if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                            continue;
                        }
                        
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    
                    // Geodesic reconstruction: constrain by mask
                    let new_val = max_neighbor.min(current_mask);
                    
                    // Use a small epsilon for floating point comparison
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }

    // ----------- Compute h-dome and create regional maxima mask (correct skimage behavior) -------------

    // Compute h-dome: this is the correct h-maxima transform
    let mut h_dome = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                h_dome[[i, j, k]] = img[[i, j, k]] - recon[[i, j, k]];
            }
        }
    }

    // Create binary mask of h-maxima regions (like skimage does)
    // These are regions where h_dome > small threshold
    let threshold = h - 1e-6; // Use h minus a small epsilon for floating point stability
    let mut h_maxima_mask = Array3::<bool>::from_elem((nx, ny, nz), false);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if img[[i, j, k]] > -1e5 && h_dome[[i, j, k]] >= threshold {
                    h_maxima_mask[[i, j, k]] = true;
                }
            }
        }
    }

    // Find peak coordinates within each h-maxima region
    // Use parallel processing with pre-allocated vectors for better performance
    let coord_vals: Vec<(Vec<[usize; 3]>, Vec<f32>)> = (0..nx)
        .into_par_iter()
        .map(|i| {
            let mut coords = Vec::new();
            let mut vals = Vec::new();

            for j in 0..ny {
                for k in 0..nz {
                    // Only process points in h-maxima regions
                    if !h_maxima_mask[[i, j, k]] {
                        continue;
                    }

                    // Check if this is a local maximum in the ORIGINAL image within h-maxima region
                    let current_val = img[[i, j, k]];
                    let mut is_local_maximum = true;
                    
                    // Check all 26 neighbors
                    'neighbor_loop: for di in -1_isize..=1_isize {
                        for dj in -1_isize..=1_isize {
                            for dk in -1_isize..=1_isize {
                                if di == 0 && dj == 0 && dk == 0 {
                                    continue; // Skip center pixel
                                }
                                
                                let ni = i as isize + di;
                                let nj = j as isize + dj;
                                let nk = k as isize + dk;

                                // Check bounds (no periodic wrapping for maxima detection)
                                if ni < 0 || nj < 0 || nk < 0 {
                                    continue;
                                }
                                let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
                                if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                                    continue;
                                }

                                // Compare ORIGINAL image values (not h-dome values)
                                // Only compare within h-maxima region OR with similar values
                                let neighbor_val = img[[ni_u, nj_u, nk_u]];
                                if neighbor_val >= current_val {
                                    is_local_maximum = false;
                                    break 'neighbor_loop;
                                }
                            }
                        }
                    }

                    if is_local_maximum {
                        coords.push([i, j, k]);
                        vals.push(current_val);
                    }
                }
            }
            (coords, vals)
        })
        .collect();

    // Combine results from all chunks with pre-allocated capacity
    let total_capacity: usize = coord_vals.iter().map(|(coords, _)| coords.len()).sum();
    let mut maxima_coords: Vec<[usize; 3]> = Vec::with_capacity(total_capacity);
    let mut maxima_vals: Vec<f32> = Vec::with_capacity(total_capacity);

    for (coords, vals) in coord_vals {
        maxima_coords.extend(coords);
        maxima_vals.extend(vals);
    }

    // Early return if no maxima found
    let n_raw = maxima_coords.len();
    if n_raw == 0 {
        let empty = Array2::<usize>::zeros((0, 3));
        return Ok((
            PyArray2::from_owned_array(py, empty.clone()).to_owned(),
            PyArray2::from_owned_array(py, empty.clone()).to_owned(),
            PyArray2::from_owned_array(py, empty).to_owned(),
        ));
    }

    // ----------- Optional peak refinement using gradient ascent -------------------
    
    if refine_peaks {
        // 27-connectivity neighbors (including center)
        let neigh_offsets: [(isize, isize, isize); 27] = [
            (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
            (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
            (-1, 1, -1), (-1, 1, 0), (-1, 1, 1),
            (0, -1, -1), (0, -1, 0), (0, -1, 1),
            (0, 0, -1), (0, 0, 0), (0, 0, 1),
            (0, 1, -1), (0, 1, 0), (0, 1, 1),
            (1, -1, -1), (1, -1, 0), (1, -1, 1),
            (1, 0, -1), (1, 0, 0), (1, 0, 1),
            (1, 1, -1), (1, 1, 0), (1, 1, 1),
        ];

        // Refine each maximum using gradient ascent
        let refined: Vec<([usize; 3], f32)> = maxima_coords
            .into_par_iter()
            .zip(maxima_vals.into_par_iter())
            .map(|(coord, _original_val)| {
                let mut pos = coord;
                let mut val = img[[pos[0], pos[1], pos[2]]];
                
                for _ in 0..max_steps {
                    let mut best = pos;
                    let mut best_val = val;
                    
                    // Check all 27 neighbors (including center for plateau detection)
                    for &(dx, dy, dz) in &neigh_offsets {
                        // Use periodic wrapping for gradient ascent
                        let ni = ((pos[0] as isize + dx).rem_euclid(nx as isize)) as usize;
                        let nj = ((pos[1] as isize + dy).rem_euclid(ny as isize)) as usize;
                        let nk = ((pos[2] as isize + dz).rem_euclid(nz as isize)) as usize;
                        let v = img[[ni, nj, nk]];
                        
                        if v > best_val {
                            best_val = v;
                            best = [ni, nj, nk];
                        }
                    }
                    
                    // Stop if we've reached a plateau or local maximum
                    if best_val <= val {
                        break;
                    }
                    pos = best;
                    val = best_val;
                }
                (pos, val)
            })
            .collect();

        // Update coordinates and values with refined results
        maxima_coords = refined.iter().map(|(coord, _)| *coord).collect();
        maxima_vals = refined.iter().map(|(_, val)| *val).collect();
    }

    // ----------- Merging with spatial indexing (periodic and non-periodic) -------------------

    let dims = (nx, ny, nz);
    let cell_ref = cell.as_ref().map(|c| c.as_array());

    let merged_pbc_coords = merge_maxima_optimized(
        &maxima_coords,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        true, // periodic wrapping for merging
    );

    let merged_cart_coords = merge_maxima_optimized(
        &maxima_coords,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        false, // no wrapping for merging
    );

    // Convert to NumPy arrays with pre-calculated sizes
    let raw_arr =
        Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords.into_iter().flatten().collect())
            .unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec(
        (merged_pbc_coords.len(), 3),
        merged_pbc_coords.into_iter().flatten().collect(),
    )
    .unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec(
        (merged_cart_coords.len(), 3),
        merged_cart_coords.into_iter().flatten().collect(),
    )
    .unwrap();

    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
}

/// Detect h-maxima after morphological opening by reconstruction.
/// This removes all maxima smaller than the structuring element (ball of radius opening_radius voxels)
/// before applying the h-maxima transform. This is more selective than Gaussian smoothing and preserves
/// the physical meaning of h-maxima. The erosion step is O(r^3) per voxel for radius r, so large radii
/// may be slow for large grids.
#[pyfunction]
#[pyo3(signature = (grid, spacing, h, opening_radius, threshold_abs=None, min_distance=0.0, cell=None, refine_peaks=false, max_steps=100))]
fn detect_maxima_h_morph_opening(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    spacing: f32,
    h: f32,
    opening_radius: usize, // radius of structuring element in grid units
    threshold_abs: Option<f32>,
    min_distance: f32,
    cell: Option<PyReadonlyArray2<'_, f64>>,
    refine_peaks: bool,
    max_steps: usize,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3};
    let mut img: Array3<f32> = grid.as_array().to_owned();
    let (nx, ny, nz) = img.dim();

    // Apply absolute threshold if provided
    if let Some(tabs) = threshold_abs {
        img.mapv_inplace(|v| if v < tabs { -1e6 } else { v });
    }

    // Debug: Print stats for true original image (before thresholding)
    let orig_data = grid.as_array();
    let mut orig_min = f32::MAX;
    let mut orig_max = f32::MIN;
    let mut orig_sum = 0.0f64;
    let mut orig_count = 0usize;
    let mut orig_above_thresh = 0usize;
    for v in orig_data.iter() {
        orig_min = orig_min.min(*v);
        orig_max = orig_max.max(*v);
        orig_sum += *v as f64;
        orig_count += 1;
        if let Some(tabs) = threshold_abs {
            if *v > tabs { orig_above_thresh += 1; }
        }
    }
    eprintln!("[DEBUG] True original image: min={:.3}, max={:.3}, mean={:.3}, count={}, above_thresh={}", orig_min, orig_max, orig_sum/(orig_count as f64), orig_count, orig_above_thresh);

    // Create mask for threshold_abs (if set)
    let mask: Option<ndarray::Array3<bool>> = threshold_abs.map(|tabs| orig_data.mapv(|v| v >= tabs));
    if let Some(ref m) = mask {
        let mask_count = m.iter().filter(|&&b| b).count();
        eprintln!("[DEBUG] Mask: {} voxels above threshold, {} below", mask_count, orig_count - mask_count);
    }

    // Copy original data for processing
    let mut img: Array3<f32> = orig_data.to_owned();
    let (nx, ny, nz) = img.dim();

    // --- Step 1: Morphological opening by reconstruction ---
    // Erode with a ball structuring element of radius opening_radius, allowing mask to shrink
    let mut eroded = img.clone();
    let mut eroded_mask = mask.clone();
    let r = opening_radius as isize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                // If masked out, skip
                if let Some(ref m) = mask {
                    if !m[[i, j, k]] {
                        if let Some(ref mut em) = eroded_mask {
                            em[[i, j, k]] = false;
                        }
                        continue;
                    }
                }
                let mut min_val = img[[i, j, k]];
                let mut fully_inside = true;
                for di in -r..=r {
                    for dj in -r..=r {
                        for dk in -r..=r {
                            if di*di + dj*dj + dk*dk > r*r { continue; }
                            let ni = i as isize + di;
                            let nj = j as isize + dj;
                            let nk = k as isize + dk;
                            if ni < 0 || nj < 0 || nk < 0 {
                                fully_inside = false;
                                break;
                            }
                            let (ni, nj, nk) = (ni as usize, nj as usize, nk as usize);
                            if ni >= nx || nj >= ny || nk >= nz {
                                fully_inside = false;
                                break;
                            }
                            if let Some(ref m) = mask {
                                if !m[[ni, nj, nk]] {
                                    fully_inside = false;
                                    break;
                                }
                            }
                            let v = img[[ni, nj, nk]];
                            if v < min_val { min_val = v; }
                        }
                        if !fully_inside { break; }
                    }
                    if !fully_inside { break; }
                }
                if fully_inside {
                    eroded[[i, j, k]] = min_val;
                    if let Some(ref mut em) = eroded_mask {
                        em[[i, j, k]] = true;
                    }
                } else {
                    if let Some(ref mut em) = eroded_mask {
                        em[[i, j, k]] = false;
                    }
                }
            }
        }
    }
    // Debug: Print stats for eroded image (unmasked voxels only)
    let mut eroded_min = f32::MAX;
    let mut eroded_max = f32::MIN;
    let mut eroded_sum = 0.0f64;
    let mut eroded_count = 0usize;
    let mut eroded_above_thresh = 0usize;
    for ((i, j, k), v) in eroded.indexed_iter() {
        if let Some(ref em) = eroded_mask {
            if !em[[i, j, k]] { continue; }
        }
        eroded_min = eroded_min.min(*v);
        eroded_max = eroded_max.max(*v);
        eroded_sum += *v as f64;
        eroded_count += 1;
        if let Some(tabs) = threshold_abs {
            if *v > tabs { eroded_above_thresh += 1; }
        }
    }
    eprintln!("[DEBUG] Eroded image: min={:.3}, max={:.3}, mean={:.3}, count={}, above_thresh={}", eroded_min, eroded_max, eroded_sum/(eroded_count as f64), eroded_count, eroded_above_thresh);
    if let Some(ref em) = eroded_mask {
        let mask_count = em.iter().filter(|&&b| b).count();
        eprintln!("[DEBUG] Eroded mask: {} voxels remain after erosion", mask_count);
    }

    // Geodesic dilation (reconstruction) of eroded under mask img, only on unmasked voxels
    let mut recon = eroded.clone();
    let mut recon_mask = eroded_mask.clone();
    let neighbors: &[(isize, isize, isize)] = &[
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
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if let Some(ref rm) = recon_mask {
                        if !rm[[i, j, k]] { continue; }
                    }
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        if ni < 0 || nj < 0 || nk < 0 { continue; }
                        let (ni, nj, nk) = (ni as usize, nj as usize, nk as usize);
                        if ni >= nx || nj >= ny || nk >= nz { continue; }
                        if let Some(ref rm) = recon_mask {
                            if !rm[[ni, nj, nk]] { continue; }
                        }
                        let v = recon[[ni, nj, nk]];
                        if v > max_neighbor { max_neighbor = v; }
                    }
                    let new_val = max_neighbor.min(img[[i, j, k]]);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    // Debug: Print stats for reconstructed (opened) image (unmasked voxels only)
    let mut recon_min = f32::MAX;
    let mut recon_max = f32::MIN;
    let mut recon_sum = 0.0f64;
    let mut recon_count = 0usize;
    let mut recon_above_thresh = 0usize;
    for ((i, j, k), v) in recon.indexed_iter() {
        if let Some(ref rm) = recon_mask {
            if !rm[[i, j, k]] { continue; }
        }
        recon_min = recon_min.min(*v);
        recon_max = recon_max.max(*v);
        recon_sum += *v as f64;
        recon_count += 1;
        if let Some(tabs) = threshold_abs {
            if *v > tabs { recon_above_thresh += 1; }
        }
    }
    eprintln!("[DEBUG] Opened image: min={:.3}, max={:.3}, mean={:.3}, count={}, above_thresh={}", recon_min, recon_max, recon_sum/(recon_count as f64), recon_count, recon_above_thresh);
    if let Some(ref rm) = recon_mask {
        let mask_count = rm.iter().filter(|&&b| b).count();
        eprintln!("[DEBUG] Opened mask: {} voxels remain after opening", mask_count);
    }

    // --- Step 2: h-maxima detection on the opened image ---
    // For clarity, call the existing detect_maxima_h_morph on recon
    let recon_py = PyArray3::from_owned_array(py, recon).to_owned();
    detect_maxima_h_morph(
        py,
        recon_py.as_ref(py).readonly(),
        spacing,
        h,
        threshold_abs,
        min_distance,
        cell,
        refine_peaks,
        max_steps,
    )
}

// -----------------------------------------------------------------------------
// Helper: merge maxima with/without periodic wrapping
// -----------------------------------------------------------------------------
fn merge_maxima(
    coords: &[[usize; 3]],
    values: &[f32],
    _spacing: f32, // Unused but kept for API compatibility
    min_distance: f32,
    cell: Option<&ndarray::ArrayView2<'_, f64>>, // full cell (Å)
    dims: (usize, usize, usize),
    periodic: bool,
) -> Vec<[usize; 3]> {
    let (nx, ny, nz) = dims;
    let n = coords.len();
    if n == 0 {
        return Vec::new();
    }

    // If min_distance is 0.0, skip merging entirely
    if min_distance <= 0.0 {
        return coords.to_vec();
    }
    let merge_radius = min_distance;

    let mut uf = UnionFind::new(n);

    for a in 0..n {
        for b in (a + 1)..n {
            let [ia, ja, ka] = coords[a];
            let [ib, jb, kb] = coords[b];

            // Convert voxel coordinates to fractional coordinates (-0.5 to 0.5)
            // This matches the coordinate system used in dgrid_from_atoms (linspace convention)
            let frac_a = [
                -0.5 + (ia as f64) / ((nx - 1) as f64),
                -0.5 + (ja as f64) / ((ny - 1) as f64),
                -0.5 + (ka as f64) / ((nz - 1) as f64),
            ];
            let frac_b = [
                -0.5 + (ib as f64) / ((nx - 1) as f64),
                -0.5 + (jb as f64) / ((ny - 1) as f64),
                -0.5 + (kb as f64) / ((nz - 1) as f64),
            ];

            // Calculate fractional difference
            let mut dfrac = [
                frac_a[0] - frac_b[0],
                frac_a[1] - frac_b[1],
                frac_a[2] - frac_b[2],
            ];

            // Apply minimum image convention in fractional space if periodic
            if periodic {
                dfrac[0] -= dfrac[0].round();
                dfrac[1] -= dfrac[1].round();
                dfrac[2] -= dfrac[2].round();
            }

            // Cartesian distance (Å) using cell matrix
            // This works for both orthogonal and non-orthogonal cells
            let dist = if let Some(cmat) = cell {
                // Transform fractional difference to Cartesian using cell matrix
                // cell matrix has lattice vectors as columns: r = cell * frac_coords
                let cart_x =
                    cmat[[0, 0]] * dfrac[0] + cmat[[0, 1]] * dfrac[1] + cmat[[0, 2]] * dfrac[2];
                let cart_y =
                    cmat[[1, 0]] * dfrac[0] + cmat[[1, 1]] * dfrac[1] + cmat[[1, 2]] * dfrac[2];
                let cart_z =
                    cmat[[2, 0]] * dfrac[0] + cmat[[2, 1]] * dfrac[1] + cmat[[2, 2]] * dfrac[2];
                ((cart_x.powi(2) + cart_y.powi(2) + cart_z.powi(2)).sqrt()) as f32
            } else {
                // No cell matrix provided - cannot compute proper distances
                // This should not happen for crystallographic systems
                return Vec::new();
            };

            if dist < merge_radius {
                uf.union(a, b);
            }
        }
    }

    // pick the member with the highest value within each cluster
    let mut merged: Vec<[usize; 3]> = Vec::new();
    let clusters = uf.get_clusters();
    for indices in clusters.values() {
        // Pick the member with highest value (no boundary preference for PBC)
        let mut best = indices[0];
        for &idx in indices {
            if values[idx] > values[best] {
                best = idx;
            }
        }
        merged.push(coords[best]);
    }
    merged
}

// -----------------------------------------------------------------------------
// Optimized merging with spatial indexing for better performance on large datasets
// -----------------------------------------------------------------------------
fn merge_maxima_optimized(
    coords: &[[usize; 3]],
    values: &[f32],
    spacing: f32, // Still used for bin_size calculation
    min_distance: f32,
    cell: Option<&ndarray::ArrayView2<'_, f64>>,
    dims: (usize, usize, usize),
    periodic: bool,
) -> Vec<[usize; 3]> {
    let (nx, ny, nz) = dims;
    let n = coords.len();
    if n == 0 {
        return Vec::new();
    }

    // If min_distance is 0.0, skip merging entirely
    if min_distance <= 0.0 {
        return coords.to_vec();
    }
    let merge_radius = min_distance;

    // For small datasets, use the simple O(n²) approach
    if n < 1000 {
        return merge_maxima(coords, values, spacing, min_distance, cell, dims, periodic);
    }

    // For larger datasets, use spatial binning to reduce comparisons
    let mut uf = UnionFind::new(n);

    // Create spatial bins to reduce pairwise comparisons
    let bin_size = (merge_radius / spacing).max(1.0) as usize;
    let bin_nx = (nx + bin_size - 1) / bin_size;
    let bin_ny = (ny + bin_size - 1) / bin_size;
    let bin_nz = (nz + bin_size - 1) / bin_size;

    // Map each coordinate to its spatial bin
    let mut bins: HashMap<(usize, usize, usize), Vec<usize>> = HashMap::new();

    for (idx, &[i, j, k]) in coords.iter().enumerate() {
        let bin_i = i / bin_size;
        let bin_j = j / bin_size;
        let bin_k = k / bin_size;
        bins.entry((bin_i, bin_j, bin_k)).or_default().push(idx);
    }

    // Only compare points within neighboring bins
    for (&(bi, bj, bk), indices) in &bins {
        // Check current bin and neighboring bins (proper 3x3x3 neighborhood)
        for di in -1_isize..=1_isize {
            for dj in -1_isize..=1_isize {
                for dk in -1_isize..=1_isize {
                    let ni = bi as isize + di;
                    let nj = bj as isize + dj;
                    let nk = bk as isize + dk;

                    // Handle periodic wrapping of bins if needed
                    let (actual_ni, actual_nj, actual_nk) = if periodic {
                        (
                            ni.rem_euclid(bin_nx as isize) as usize,
                            nj.rem_euclid(bin_ny as isize) as usize,
                            nk.rem_euclid(bin_nz as isize) as usize,
                        )
                    } else {
                        if ni < 0 || nj < 0 || nk < 0 || ni >= bin_nx as isize || nj >= bin_ny as isize || nk >= bin_nz as isize {
                            continue;
                        }
                        (ni as usize, nj as usize, nk as usize)
                    };

                    if let Some(neighbor_indices) = bins.get(&(actual_ni, actual_nj, actual_nk)) {
                        // Compare all pairs between current bin and neighbor bin
                        for &a in indices {
                            for &b in neighbor_indices {
                                if a >= b {
                                    continue;
                                } // Avoid duplicate comparisons

                                let [ia, ja, ka] = coords[a];
                                let [ib, jb, kb] = coords[b];

                                // Convert voxel coordinates to fractional coordinates (-0.5 to 0.5)
                                // This matches the coordinate system used in dgrid_from_atoms (linspace convention)
                                let frac_a = [
                                    -0.5 + (ia as f64) / ((nx - 1) as f64),
                                    -0.5 + (ja as f64) / ((ny - 1) as f64),
                                    -0.5 + (ka as f64) / ((nz - 1) as f64),
                                ];
                                let frac_b = [
                                    -0.5 + (ib as f64) / ((nx - 1) as f64),
                                    -0.5 + (jb as f64) / ((ny - 1) as f64),
                                    -0.5 + (kb as f64) / ((nz - 1) as f64),
                                ];

                                // Calculate fractional difference
                                let mut dfrac = [
                                    frac_a[0] - frac_b[0],
                                    frac_a[1] - frac_b[1],
                                    frac_a[2] - frac_b[2],
                                ];

                                // Apply minimum image convention in fractional space if periodic
                                if periodic {
                                    dfrac[0] -= dfrac[0].round();
                                    dfrac[1] -= dfrac[1].round();
                                    dfrac[2] -= dfrac[2].round();
                                }

                                // Cartesian distance using cell matrix
                                // This works for both orthogonal and non-orthogonal cells
                                let dist = if let Some(cmat) = cell {
                                    // Transform fractional difference to Cartesian using cell matrix
                                    // cell matrix has lattice vectors as columns: r = cell * frac_coords
                                    let cart_x = cmat[[0, 0]] * dfrac[0]
                                        + cmat[[0, 1]] * dfrac[1]
                                        + cmat[[0, 2]] * dfrac[2];
                                    let cart_y = cmat[[1, 0]] * dfrac[0]
                                        + cmat[[1, 1]] * dfrac[1]
                                        + cmat[[1, 2]] * dfrac[2];
                                    let cart_z = cmat[[2, 0]] * dfrac[0]
                                        + cmat[[2, 1]] * dfrac[1]
                                        + cmat[[2, 2]] * dfrac[2];
                                    ((cart_x.powi(2) + cart_y.powi(2) + cart_z.powi(2)).sqrt())
                                        as f32
                                } else {
                                    // No cell matrix provided - cannot compute proper distances
                                    // Skip this comparison
                                    continue;
                                };

                                if dist < merge_radius {
                                    uf.union(a, b);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Extract clusters and pick best representative from each
    let mut merged: Vec<[usize; 3]> = Vec::new();
    let clusters = uf.get_clusters();
    for indices in clusters.values() {
        // Pick the member with highest value (no boundary preference for PBC)
        let mut best = indices[0];
        for &idx in indices {
            if values[idx] > values[best] {
                best = idx;
            }
        }
        merged.push(coords[best]);
    }
    merged
}

/// Move maxima to exact peak locations using gradient ascent
/// This refines the peak locations found by h-maxima detection
#[pyfunction]
#[pyo3(signature = (grid, starts, max_steps=100))]
fn move_maxima_to_peaks(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    starts: PyReadonlyArray2<'_, usize>, // (n,3)
    max_steps: Option<usize>,
) -> PyResult<(Py<PyArray2<usize>>, Py<PyArray1<f32>>)> {
    use ndarray::{Array1, Array2};
    use rayon::prelude::*;
    
    // Clone into owned ndarrays so the underlying data is independent of the
    // PyArray and therefore `Sync`/`Send`.
    let g = grid.as_array().to_owned();
    let starts_owned = starts.as_array().to_owned();

    let (nx, ny, nz) = g.dim();
    let n = starts_owned.dim().0;
    let ms = max_steps.unwrap_or(100);

    // 27-connectivity neighbors (including center)
    let neigh_offsets: [(isize, isize, isize); 27] = [
        (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
        (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1),
        (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1), (0, 0, 0), (0, 0, 1),
        (0, 1, -1), (0, 1, 0), (0, 1, 1),
        (1, -1, -1), (1, -1, 0), (1, -1, 1),
        (1, 0, -1), (1, 0, 0), (1, 0, 1),
        (1, 1, -1), (1, 1, 0), (1, 1, 1),
    ];

    let peaks: Vec<([usize; 3], f32)> = (0..n)
        .into_par_iter()
        .map(|idx| {
            let s_row = starts_owned.row(idx);
            let mut pos = [s_row[0], s_row[1], s_row[2]];
            let mut val = g[[pos[0], pos[1], pos[2]]];
            
            for _ in 0..ms {
                let mut best = pos;
                let mut best_val = val;
                
                // Check all 27 neighbors (including center for plateau detection)
                for &(dx, dy, dz) in &neigh_offsets {
                    // Use periodic wrapping for gradient ascent
                    let ni = ((pos[0] as isize + dx).rem_euclid(nx as isize)) as usize;
                    let nj = ((pos[1] as isize + dy).rem_euclid(ny as isize)) as usize;
                    let nk = ((pos[2] as isize + dz).rem_euclid(nz as isize)) as usize;
                    let v = g[[ni, nj, nk]];
                    
                    if v > best_val {
                        best_val = v;
                        best = [ni, nj, nk];
                    }
                }
                
                // Stop if we've reached a plateau or local maximum
                if best_val <= val {
                    break;
                }
                pos = best;
                val = best_val;
            }
            (pos, val)
        })
        .collect();

    // Split into separate arrays
    let mut coords_vec: Vec<usize> = Vec::with_capacity(n * 3);
    let mut vals_vec: Vec<f32> = Vec::with_capacity(n);
    for (p, v) in peaks {
        coords_vec.extend_from_slice(&p);
        vals_vec.push(v);
    }
    let coords_arr = Array2::from_shape_vec((n, 3), coords_vec).expect("shape error");
    let vals_arr = Array1::from_vec(vals_vec);
    
    Ok((
        PyArray2::from_owned_array(py, coords_arr).to_owned(),
        PyArray1::from_owned_array(py, vals_arr).to_owned(),
    ))
}

/// Detect h-maxima with separate mask_thickness and threshold_abs, and configurable neighborhood size.
/// Uses periodic boundary conditions throughout the entire detection process.
/// - mask_thickness: mask out voxels with value < mask_thickness before maxima detection.
/// - threshold_abs: post-filter maxima by value after detection.
/// - neighborhood_size: size of the cube neighborhood for local maxima (default 5).
#[pyfunction]
#[pyo3(signature = (grid, spacing, h, mask_thickness=None, threshold_abs=None, min_distance=0.0, cell=None, refine_peaks=false, max_steps=100, neighborhood_size=5))]
fn detect_maxima_h_morph2(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    spacing: f32,
    h: f32,
    mask_thickness: Option<f32>,
    threshold_abs: Option<f32>,
    min_distance: f32,
    cell: Option<PyReadonlyArray2<'_, f64>>,
    refine_peaks: bool,
    max_steps: usize,
    neighborhood_size: usize,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3, Array2};
    use rayon::prelude::*;
    let mut img: Array3<f32> = grid.as_array().to_owned();
    let (nx, ny, nz) = img.dim();

    // Step 1: Mask out voxels with value < mask_thickness (if set)
    let mask: Option<Array3<bool>> = mask_thickness.map(|mt| img.mapv(|v| v >= mt));

    // Step 2: h-maxima transform with periodic boundary conditions
    if let Some(mt) = mask_thickness {
        img.mapv_inplace(|v| if v < mt { -1e6 } else { v });
    }
    // Create marker image (f - h)
    let marker = img.mapv(|v| v - h);
    let mut recon = marker.clone();
    // 26-connectivity neighbors for reconstruction (with periodic wrapping)
    let neighbors: &[(isize, isize, isize)] = &[
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
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if let Some(ref m) = mask {
                        if !m[[i, j, k]] { continue; }
                    }
                    let current_mask = img[[i, j, k]];
                    if current_mask <= -1e5 { continue; }
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        // Apply periodic boundary conditions using rem_euclid
                        let ni_u = ni.rem_euclid(nx as isize) as usize;
                        let nj_u = nj.rem_euclid(ny as isize) as usize;
                        let nk_u = nk.rem_euclid(nz as isize) as usize;
                        if let Some(ref m) = mask {
                            if !m[[ni_u, nj_u, nk_u]] { continue; }
                        }
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    let new_val = max_neighbor.min(current_mask);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    // Compute h-dome
    let mut h_dome = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                h_dome[[i, j, k]] = img[[i, j, k]] - recon[[i, j, k]];
            }
        }
    }
    // h-maxima mask
    let threshold = h - 1e-6;
    let mut h_maxima_mask = Array3::<bool>::from_elem((nx, ny, nz), false);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if img[[i, j, k]] > -1e5 && h_dome[[i, j, k]] >= threshold {
                    h_maxima_mask[[i, j, k]] = true;
                }
            }
        }
    }
    // Step 3: Find local maxima in a larger neighborhood with periodic boundary conditions
    let nsize = neighborhood_size as isize;
    let half = nsize / 2;
    let mut maxima_coords = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if !h_maxima_mask[[i, j, k]] { continue; }
                let current_val = img[[i, j, k]];
                let mut is_local_maximum = true;
                'neighbor_loop: for di in -half..=half {
                    for dj in -half..=half {
                        for dk in -half..=half {
                            if di == 0 && dj == 0 && dk == 0 { continue; }
                            let ni = i as isize + di;
                            let nj = j as isize + dj;
                            let nk = k as isize + dk;
                            // Apply periodic boundary conditions using rem_euclid
                            let ni_u = ni.rem_euclid(nx as isize) as usize;
                            let nj_u = nj.rem_euclid(ny as isize) as usize;
                            let nk_u = nk.rem_euclid(nz as isize) as usize;
                            if let Some(ref m) = mask {
                                if !m[[ni_u, nj_u, nk_u]] { continue; }
                            }
                            let neighbor_val = img[[ni_u, nj_u, nk_u]];
                            if neighbor_val >= current_val {
                                is_local_maximum = false;
                                break 'neighbor_loop;
                            }
                        }
                    }
                }
                if is_local_maximum {
                    maxima_coords.push([i, j, k]);
                }
            }
        }
    }
    // Step 4: Post-filter maxima by threshold_abs (if set)
    let mut maxima_coords_filtered = Vec::new();
    for coord in &maxima_coords {
        let val = img[[coord[0], coord[1], coord[2]]];
        if let Some(tabs) = threshold_abs {
            if val < tabs { continue; }
        }
        maxima_coords_filtered.push(*coord);
    }
    // Step 5: Merge maxima (periodic and non-periodic)
    let cell_ref = cell.as_ref().map(|c| c.as_array());
    let dims = (nx, ny, nz);
    let maxima_vals: Vec<f32> = maxima_coords_filtered.iter().map(|c| img[[c[0], c[1], c[2]]]).collect();
    let merged_pbc_coords = merge_maxima_optimized(
        &maxima_coords_filtered,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        true,
    );
    let merged_cart_coords = merge_maxima_optimized(
        &maxima_coords_filtered,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        false,
    );
    // Convert to NumPy arrays
    let n_raw = maxima_coords_filtered.len();
    let raw_arr = Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords_filtered.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec((merged_pbc_coords.len(), 3), merged_pbc_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec((merged_cart_coords.len(), 3), merged_cart_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
}

/// Simplified h-maxima detection with periodic boundary conditions.
/// Removes the strict local maxima filtering for better maxima detection.
/// - mask_thickness: mask out voxels with value < mask_thickness before maxima detection.
/// - threshold_abs: post-filter maxima by value after detection.
/// - neighborhood_size: kept for API compatibility (not used in simplified version).
/// - h_dome_tolerance: tolerance for h-dome threshold (default 0.1*h, looser than original).
#[pyfunction]
#[pyo3(signature = (grid, spacing, h, mask_thickness=None, threshold_abs=None, min_distance=0.0, cell=None, refine_peaks=false, max_steps=100, neighborhood_size=5, h_dome_tolerance=None))]
fn detect_maxima_h_morph_simple(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    spacing: f32,
    h: f32,
    mask_thickness: Option<f32>,
    threshold_abs: Option<f32>,
    min_distance: f32,
    cell: Option<PyReadonlyArray2<'_, f64>>,
    refine_peaks: bool,
    max_steps: usize,
    neighborhood_size: usize, // Kept for API compatibility
    h_dome_tolerance: Option<f32>,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3, Array2};
    let mut img: Array3<f32> = grid.as_array().to_owned();
    let (nx, ny, nz) = img.dim();

    // Step 1: Mask out voxels with value < mask_thickness (if set)
    let mask: Option<Array3<bool>> = mask_thickness.map(|mt| img.mapv(|v| v >= mt));

    // Step 2: h-maxima transform with periodic boundary conditions
    if let Some(mt) = mask_thickness {
        img.mapv_inplace(|v| if v < mt { -1e6 } else { v });
    }
    // Create marker image (f - h)
    let marker = img.mapv(|v| v - h);
    let mut recon = marker.clone();
    // 26-connectivity neighbors for reconstruction (with periodic wrapping)
    let neighbors: &[(isize, isize, isize)] = &[
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
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if let Some(ref m) = mask {
                        if !m[[i, j, k]] { continue; }
                    }
                    let current_mask = img[[i, j, k]];
                    if current_mask <= -1e5 { continue; }
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        // Apply periodic boundary conditions using rem_euclid
                        let ni_u = ni.rem_euclid(nx as isize) as usize;
                        let nj_u = nj.rem_euclid(ny as isize) as usize;
                        let nk_u = nk.rem_euclid(nz as isize) as usize;
                        if let Some(ref m) = mask {
                            if !m[[ni_u, nj_u, nk_u]] { continue; }
                        }
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    let new_val = max_neighbor.min(current_mask);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    
    // Step 3: Compute h-dome with looser threshold
    let mut h_dome = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                h_dome[[i, j, k]] = img[[i, j, k]] - recon[[i, j, k]];
            }
        }
    }
    
    // Use looser h-dome threshold (default 0.1*h instead of h - 1e-6)
    let tolerance = h_dome_tolerance.unwrap_or(0.1 * h);
    let threshold = h - tolerance;
    let mut maxima_coords = Vec::new();
    
    // Step 4 REMOVED: Collect all h-maxima points without additional local maxima filtering
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if let Some(ref m) = mask {
                    if !m[[i, j, k]] { continue; }
                }
                if img[[i, j, k]] > -1e5 && h_dome[[i, j, k]] >= threshold {
                    maxima_coords.push([i, j, k]);
                }
            }
        }
    }

    // Step 5: Post-filter maxima by threshold_abs (if set)
    let mut maxima_coords_filtered = Vec::new();
    for coord in &maxima_coords {
        let val = img[[coord[0], coord[1], coord[2]]];
        if let Some(tabs) = threshold_abs {
            if val < tabs { continue; }
        }
        maxima_coords_filtered.push(*coord);
    }
    
    // Step 6: Merge maxima (periodic and non-periodic)
    let cell_ref = cell.as_ref().map(|c| c.as_array());
    let dims = (nx, ny, nz);
    let maxima_vals: Vec<f32> = maxima_coords_filtered.iter().map(|c| img[[c[0], c[1], c[2]]]).collect();
    let merged_pbc_coords = merge_maxima_optimized(
        &maxima_coords_filtered,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        true,
    );
    let merged_cart_coords = merge_maxima_optimized(
        &maxima_coords_filtered,
        &maxima_vals,
        spacing,
        min_distance,
        cell_ref.as_ref(),
        dims,
        false,
    );
    
    // Convert to NumPy arrays
    let n_raw = maxima_coords_filtered.len();
    let raw_arr = Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords_filtered.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec((merged_pbc_coords.len(), 3), merged_pbc_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec((merged_cart_coords.len(), 3), merged_cart_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    
    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
}

/// Simplified h-maxima detection with periodic boundary conditions.
/// Removes spacing parameter - uses direct pairwise PBC distance calculations.
/// - mask_thickness: mask out voxels with value < mask_thickness before maxima detection.
/// - threshold_abs: post-filter maxima by value after detection.
/// - h_dome_tolerance: tolerance for h-dome threshold (default 0.1*h, looser than original).
#[pyfunction]
#[pyo3(signature = (grid, h, mask_thickness=None, threshold_abs=None, min_distance=0.0, cell=None, h_dome_tolerance=None))]
fn detect_maxima_h_morph_no_spacing(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    h: f32,
    mask_thickness: Option<f32>,
    threshold_abs: Option<f32>,
    min_distance: f32,
    cell: Option<PyReadonlyArray2<'_, f64>>,
    h_dome_tolerance: Option<f32>,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3, Array2};
    let mut img: Array3<f32> = grid.as_array().to_owned();
    let (nx, ny, nz) = img.dim();

    // Step 1: Mask out voxels with value < mask_thickness (if set)
    let mask: Option<Array3<bool>> = mask_thickness.map(|mt| img.mapv(|v| v >= mt));

    // Step 2: h-maxima transform with periodic boundary conditions
    if let Some(mt) = mask_thickness {
        img.mapv_inplace(|v| if v < mt { -1e6 } else { v });
    }
    
    // Create marker image (f - h)
    let marker = img.mapv(|v| v - h);
    let mut recon = marker.clone();
    
    // 26-connectivity neighbors for reconstruction (with periodic wrapping)
    let neighbors: &[(isize, isize, isize)] = &[
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
    
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    if let Some(ref m) = mask {
                        if !m[[i, j, k]] { continue; }
                    }
                    let current_mask = img[[i, j, k]];
                    if current_mask <= -1e5 { continue; }
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        // Apply periodic boundary conditions using rem_euclid
                        let ni_u = ni.rem_euclid(nx as isize) as usize;
                        let nj_u = nj.rem_euclid(ny as isize) as usize;
                        let nk_u = nk.rem_euclid(nz as isize) as usize;
                        if let Some(ref m) = mask {
                            if !m[[ni_u, nj_u, nk_u]] { continue; }
                        }
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    let new_val = max_neighbor.min(current_mask);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    
    // Step 3: Compute h-dome with looser threshold
    let mut h_dome = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                h_dome[[i, j, k]] = img[[i, j, k]] - recon[[i, j, k]];
            }
        }
    }
    
    // Use looser h-dome threshold (default 0.1*h instead of h - 1e-6)
    let tolerance = h_dome_tolerance.unwrap_or(0.1 * h);
    let threshold = h - tolerance;
    let mut maxima_coords = Vec::new();
    
    // Collect all h-maxima points without additional local maxima filtering
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if let Some(ref m) = mask {
                    if !m[[i, j, k]] { continue; }
                }
                if img[[i, j, k]] > -1e5 && h_dome[[i, j, k]] >= threshold {
                    maxima_coords.push([i, j, k]);
                }
            }
        }
    }

    // Step 4: Post-filter maxima by threshold_abs (if set)
    let mut maxima_coords_filtered = Vec::new();
    for coord in &maxima_coords {
        let val = img[[coord[0], coord[1], coord[2]]];
        if let Some(tabs) = threshold_abs {
            if val < tabs { continue; }
        }
        maxima_coords_filtered.push(*coord);
    }
    
    // Step 5: Simple pairwise distance-based merging (no spacing needed)
    let cell_ref = cell.as_ref().map(|c| c.as_array());
    let maxima_vals: Vec<f32> = maxima_coords_filtered.iter().map(|c| img[[c[0], c[1], c[2]]]).collect();
    
    let merged_pbc_coords = merge_maxima_simple(
        &maxima_coords_filtered,
        &maxima_vals,
        min_distance,
        cell_ref.as_ref(),
        (nx, ny, nz),
        true, // periodic
    );
    
    let merged_cart_coords = merge_maxima_simple(
        &maxima_coords_filtered,
        &maxima_vals,
        min_distance,
        cell_ref.as_ref(),
        (nx, ny, nz),
        false, // non-periodic
    );
    
    // Convert to NumPy arrays
    let n_raw = maxima_coords_filtered.len();
    let raw_arr = Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords_filtered.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec((merged_pbc_coords.len(), 3), merged_pbc_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec((merged_cart_coords.len(), 3), merged_cart_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    
    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
}

/// Simple merging function without spacing parameter
fn merge_maxima_simple(
    coords: &[[usize; 3]],
    values: &[f32],
    min_distance: f32,
    cell: Option<&ndarray::ArrayView2<'_, f64>>,
    dims: (usize, usize, usize),
    periodic: bool,
) -> Vec<[usize; 3]> {
    let (nx, ny, nz) = dims;
    let n = coords.len();
    if n == 0 {
        return Vec::new();
    }

    // If min_distance is 0.0, skip merging entirely
    if min_distance <= 0.0 {
        return coords.to_vec();
    }
    let merge_radius = min_distance;

    let mut uf = UnionFind::new(n);

    for a in 0..n {
        for b in (a + 1)..n {
            let [ia, ja, ka] = coords[a];
            let [ib, jb, kb] = coords[b];

            // Convert voxel coordinates to fractional coordinates (-0.5 to 0.5)
            let frac_a = [
                -0.5 + (ia as f64) / ((nx - 1) as f64),
                -0.5 + (ja as f64) / ((ny - 1) as f64),
                -0.5 + (ka as f64) / ((nz - 1) as f64),
            ];
            let frac_b = [
                -0.5 + (ib as f64) / ((nx - 1) as f64),
                -0.5 + (jb as f64) / ((ny - 1) as f64),
                -0.5 + (kb as f64) / ((nz - 1) as f64),
            ];

            // Calculate fractional difference
            let mut dfrac = [
                frac_a[0] - frac_b[0],
                frac_a[1] - frac_b[1],
                frac_a[2] - frac_b[2],
            ];

            // Apply minimum image convention in fractional space if periodic
            if periodic {
                dfrac[0] -= dfrac[0].round();
                dfrac[1] -= dfrac[1].round();
                dfrac[2] -= dfrac[2].round();
            }

            // Cartesian distance using cell matrix (preferred) or grid units (fallback)
            let dist = if let Some(cmat) = cell {
                // Transform fractional difference to Cartesian using cell matrix
                let cart_x = cmat[[0, 0]] * dfrac[0] + cmat[[0, 1]] * dfrac[1] + cmat[[0, 2]] * dfrac[2];
                let cart_y = cmat[[1, 0]] * dfrac[0] + cmat[[1, 1]] * dfrac[1] + cmat[[1, 2]] * dfrac[2];
                let cart_z = cmat[[2, 0]] * dfrac[0] + cmat[[2, 1]] * dfrac[1] + cmat[[2, 2]] * dfrac[2];
                ((cart_x.powi(2) + cart_y.powi(2) + cart_z.powi(2)).sqrt()) as f32
            } else {
                // Fallback: use grid units (assumes cubic grid, less accurate)
                let grid_dist = (dfrac[0].powi(2) + dfrac[1].powi(2) + dfrac[2].powi(2)).sqrt() as f32;
                grid_dist * ((nx + ny + nz) as f32 / 3.0) // rough scaling
            };

            if dist < merge_radius {
                uf.union(a, b);
            }
        }
    }

    // Pick the member with the highest value within each cluster
    let mut merged: Vec<[usize; 3]> = Vec::new();
    let clusters = uf.get_clusters();
    for indices in clusters.values() {
        let mut best = indices[0];
        for &idx in indices {
            if values[idx] > values[best] {
                best = idx;
            }
        }
        merged.push(coords[best]);
    }
    merged
}

/// Geometry-adaptive h-maxima detection using distance field properties.
/// Uses local dgrid values and gradients to adaptively choose h values:
/// - Large smooth regions (channels) use large h to suppress surface roughness
/// - Small regions (pockets) use small h for sensitive detection  
/// - High gradient regions (near surfaces) use very small h
/// 
/// Parameters:
/// - h_small: h value for small pockets and surface regions (default 0.2)
/// - h_large: h value for large channel interiors (default 0.8)  
/// - transition_radius: dgrid value separating small/large regions in Å (default 3.0)
/// - gradient_threshold: gradient magnitude indicating surface proximity (default 0.5)
#[pyfunction]
#[pyo3(signature = (grid, h_small=0.2, h_large=0.8, transition_radius=3.0, gradient_threshold=0.5, cell=None, min_distance=0.0, mask_thickness=None, threshold_abs=None))]
fn detect_maxima_adaptive_geometry(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    h_small: f32,      // For pockets (default 0.2)
    h_large: f32,      // For channels (default 0.8) 
    transition_radius: f32, // Transition point in Å (default 3.0)
    gradient_threshold: f32, // High gradient = surface (default 0.5)
    cell: Option<PyReadonlyArray2<'_, f64>>,
    min_distance: f32,
    mask_thickness: Option<f32>,
    threshold_abs: Option<f32>,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>, 
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3, Array2};
    let dgrid = grid.as_array();
    let (nx, ny, nz) = dgrid.dim();
    
    // Step 1: Compute adaptive h field based on local geometry
    let mut adaptive_h = Array3::<f32>::zeros((nx, ny, nz));
    
    for i in 1..nx-1 {
        for j in 1..ny-1 {
            for k in 1..nz-1 {
                let local_radius = dgrid[[i, j, k]];
                
                // Compute local gradient magnitude (indicates proximity to surface)
                let grad_x = (dgrid[[i+1, j, k]] - dgrid[[i-1, j, k]]) / 2.0;
                let grad_y = (dgrid[[i, j+1, k]] - dgrid[[i, j-1, k]]) / 2.0;
                let grad_z = (dgrid[[i, j, k+1]] - dgrid[[i, j, k-1]]) / 2.0;
                let gradient_mag = (grad_x*grad_x + grad_y*grad_y + grad_z*grad_z).sqrt();
                
                // Adaptive h based on local geometry
                let h = if gradient_mag > gradient_threshold {
                    // Near surface - use very small h for sensitivity to small features
                    h_small * 0.5
                } else if local_radius > transition_radius {
                    // Large channel interior - use large h to suppress surface roughness
                    h_large
                } else {
                    // Small pocket interior - use small h for detection
                    h_small
                };
                
                adaptive_h[[i, j, k]] = h;
            }
        }
    }
    
    // Handle boundaries by copying nearest interior values
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1 {
                    let ii = if i == 0 { 1 } else if i == nx-1 { nx-2 } else { i };
                    let jj = if j == 0 { 1 } else if j == ny-1 { ny-2 } else { j };
                    let kk = if k == 0 { 1 } else if k == nz-1 { nz-2 } else { k };
                    adaptive_h[[i, j, k]] = adaptive_h[[ii, jj, kk]];
                }
            }
        }
    }
    
    // Step 2: Apply adaptive h-maxima transform
    let mut img = dgrid.to_owned();
    if let Some(mt) = mask_thickness {
        img.mapv_inplace(|v| if v < mt { -1e6 } else { v });
    }
    
    // Create spatially-varying marker image: marker[i,j,k] = img[i,j,k] - adaptive_h[i,j,k]
    let mut marker = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                marker[[i, j, k]] = img[[i, j, k]] - adaptive_h[[i, j, k]];
            }
        }
    }
    
    // Step 3: Morphological reconstruction with periodic boundaries
    let mut recon = marker.clone();
    let neighbors: &[(isize, isize, isize)] = &[
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
    
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let current_mask = img[[i, j, k]];
                    if current_mask <= -1e5 { continue; }
                    
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        // Periodic boundary conditions
                        let ni_u = ni.rem_euclid(nx as isize) as usize;
                        let nj_u = nj.rem_euclid(ny as isize) as usize;
                        let nk_u = nk.rem_euclid(nz as isize) as usize;
                        
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    
                    let new_val = max_neighbor.min(current_mask);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    
    // Step 4: Compute adaptive h-dome and find maxima
    let mut maxima_coords = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if img[[i, j, k]] <= -1e5 { continue; }
                
                let h_dome_val = img[[i, j, k]] - recon[[i, j, k]];
                let local_h = adaptive_h[[i, j, k]];
                let threshold = local_h * 0.9; // Slightly looser than exact h
                
                if h_dome_val >= threshold {
                    maxima_coords.push([i, j, k]);
                }
            }
        }
    }
    
    // Step 5: Apply remaining filters and merge
    let mut maxima_coords_filtered = Vec::new();
    for coord in &maxima_coords {
        let val = img[[coord[0], coord[1], coord[2]]];
        if let Some(tabs) = threshold_abs {
            if val < tabs { continue; }
        }
        maxima_coords_filtered.push(*coord);
    }
    
    // Step 6: Merge with PBC-aware distances
    let cell_ref = cell.as_ref().map(|c| c.as_array());
    let maxima_vals: Vec<f32> = maxima_coords_filtered.iter().map(|c| img[[c[0], c[1], c[2]]]).collect();
    
    let merged_pbc_coords = merge_maxima_simple(
        &maxima_coords_filtered,
        &maxima_vals,
        min_distance,
        cell_ref.as_ref(),
        (nx, ny, nz),
        true,
    );
    
    let merged_cart_coords = merge_maxima_simple(
        &maxima_coords_filtered,
        &maxima_vals,
        min_distance,
        cell_ref.as_ref(),
        (nx, ny, nz),
        false,
    );
    
    // Convert to NumPy arrays
    let n_raw = maxima_coords_filtered.len();
    let raw_arr = Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords_filtered.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec((merged_pbc_coords.len(), 3), merged_pbc_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec((merged_cart_coords.len(), 3), merged_cart_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    
    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
}

/// Selective smoothing h-maxima detection.
/// Applies spatially-varying smoothing based on local pore size:
/// - Small pockets (dgrid < preserve_threshold): minimal smoothing to preserve features
/// - Large channels (dgrid > smooth_threshold): heavy smoothing to remove surface roughness
/// - Intermediate regions: gradual transition
/// 
/// Parameters:
/// - preserve_threshold: dgrid value below which features are preserved (default 4.0 Å)
/// - smooth_threshold: dgrid value above which heavy smoothing is applied (default 6.0 Å)  
/// - max_smooth_radius: maximum smoothing radius in grid units (default 2.0)
/// - h: single h parameter for maxima detection after smoothing (default 0.2)
/// - cell: REQUIRED for merging - unit cell matrix for distance calculations
/// - min_distance: merge distance in Å (default 1.0) - set to 0.0 to disable merging
#[pyfunction]
#[pyo3(signature = (grid, preserve_threshold=4.0, smooth_threshold=6.0, max_smooth_radius=2.0, h=0.2, cell=None, min_distance=1.0, mask_thickness=None, threshold_abs=None))]
fn detect_maxima_selective_smoothing(
    py: Python<'_>,
    grid: PyReadonlyArray3<'_, f32>,
    preserve_threshold: f32,   // Preserve features smaller than this (default 4.0 Å)
    smooth_threshold: f32,     // Smooth features larger than this (default 6.0 Å)
    max_smooth_radius: f32,    // Maximum smoothing radius in grid units (default 2.0)
    h: f32,                    // Single h parameter after smoothing (default 0.2)
    cell: Option<PyReadonlyArray2<'_, f64>>,
    min_distance: f32,
    mask_thickness: Option<f32>,
    threshold_abs: Option<f32>,
) -> PyResult<(
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
    Py<PyArray2<usize>>,
)> {
    use ndarray::{Array3, Array2};
    use std::f32::consts::PI;
    
    let dgrid = grid.as_array();
    let (nx, ny, nz) = dgrid.dim();
    
    // Step 1: Apply selective smoothing
    let mut smoothed_grid = dgrid.to_owned();
    
    // Create smoothing radius field based on local dgrid values
    let mut smooth_radii = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let local_pore_size = dgrid[[i, j, k]];
                
                let smooth_radius = if local_pore_size < preserve_threshold {
                    // Small pockets - minimal smoothing (preserve features)
                    0.0
                } else if local_pore_size > smooth_threshold {
                    // Large channels - maximum smoothing (remove roughness)
                    max_smooth_radius
                } else {
                    // Intermediate - linear interpolation
                    let fraction = (local_pore_size - preserve_threshold) / (smooth_threshold - preserve_threshold);
                    fraction * max_smooth_radius
                };
                
                smooth_radii[[i, j, k]] = smooth_radius;
            }
        }
    }
    
    // Apply spatially-varying Gaussian smoothing
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let smooth_radius = smooth_radii[[i, j, k]];
                
                if smooth_radius > 0.1 { // Only smooth if radius is significant
                    let kernel_size = (smooth_radius * 3.0).ceil() as usize; // 3-sigma kernel
                    let mut weighted_sum = 0.0f32;
                    let mut weight_total = 0.0f32;
                    
                    // Apply Gaussian kernel with periodic boundaries
                    for di in -(kernel_size as isize)..=(kernel_size as isize) {
                        for dj in -(kernel_size as isize)..=(kernel_size as isize) {
                            for dk in -(kernel_size as isize)..=(kernel_size as isize) {
                                let ni = (i as isize + di).rem_euclid(nx as isize) as usize;
                                let nj = (j as isize + dj).rem_euclid(ny as isize) as usize;
                                let nk = (k as isize + dk).rem_euclid(nz as isize) as usize;
                                
                                let dist_sq = (di * di + dj * dj + dk * dk) as f32;
                                let sigma_sq = smooth_radius * smooth_radius;
                                
                                if sigma_sq > 0.0 {
                                    let weight = (-0.5 * dist_sq / sigma_sq).exp();
                                    weighted_sum += weight * dgrid[[ni, nj, nk]];
                                    weight_total += weight;
                                }
                            }
                        }
                    }
                    
                    if weight_total > 0.0 {
                        smoothed_grid[[i, j, k]] = weighted_sum / weight_total;
                    }
                }
            }
        }
    }
    
    // Step 2: Apply standard h-maxima detection on smoothed grid
    let mut img = smoothed_grid;
    if let Some(mt) = mask_thickness {
        img.mapv_inplace(|v| if v < mt { -1e6 } else { v });
    }
    
    // Create marker image (f - h) - single h value now
    let marker = img.mapv(|v| v - h);
    let mut recon = marker.clone();
    
    // 26-connectivity neighbors for reconstruction (with periodic wrapping)
    let neighbors: &[(isize, isize, isize)] = &[
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
    
    let mut changed = true;
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;
    while changed && iteration < MAX_ITERATIONS {
        changed = false;
        iteration += 1;
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let current_mask = img[[i, j, k]];
                    if current_mask <= -1e5 { continue; }
                    
                    let mut max_neighbor = recon[[i, j, k]];
                    for &(di, dj, dk) in neighbors {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        let nk = k as isize + dk;
                        // Periodic boundary conditions
                        let ni_u = ni.rem_euclid(nx as isize) as usize;
                        let nj_u = nj.rem_euclid(ny as isize) as usize;
                        let nk_u = nk.rem_euclid(nz as isize) as usize;
                        
                        let neighbor_recon = recon[[ni_u, nj_u, nk_u]];
                        if neighbor_recon > max_neighbor {
                            max_neighbor = neighbor_recon;
                        }
                    }
                    
                    let new_val = max_neighbor.min(current_mask);
                    if (new_val - recon[[i, j, k]]).abs() > 1e-6 {
                        recon[[i, j, k]] = new_val;
                        changed = true;
                    }
                }
            }
        }
    }
    
    // Step 3: Compute h-dome and find maxima (standard approach)
    let mut h_dome = Array3::<f32>::zeros((nx, ny, nz));
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                h_dome[[i, j, k]] = img[[i, j, k]] - recon[[i, j, k]];
            }
        }
    }
    
    let threshold = h * 0.9; // Slightly looser threshold
    let mut maxima_coords = Vec::new();
    
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if img[[i, j, k]] > -1e5 && h_dome[[i, j, k]] >= threshold {
                    maxima_coords.push([i, j, k]);
                }
            }
        }
    }
    
    // Step 4: Apply remaining filters and merge
    let mut maxima_coords_filtered = Vec::new();
    for coord in &maxima_coords {
        let val = img[[coord[0], coord[1], coord[2]]];
        if let Some(tabs) = threshold_abs {
            if val < tabs { continue; }
        }
        maxima_coords_filtered.push(*coord);
    }
    
    // Step 5: Merge with PBC-aware distances  
    let cell_ref = cell.as_ref().map(|c| c.as_array());
    let maxima_vals: Vec<f32> = maxima_coords_filtered.iter().map(|c| img[[c[0], c[1], c[2]]]).collect();
    
    // Apply merging only if min_distance > 0
    let merged_pbc_coords = if min_distance > 0.0 {
        merge_maxima_simple(
            &maxima_coords_filtered,
            &maxima_vals,
            min_distance,
            cell_ref.as_ref(),
            (nx, ny, nz),
            true, // periodic
        )
    } else {
        maxima_coords_filtered.clone() // No merging
    };
    
    let merged_cart_coords = if min_distance > 0.0 {
        merge_maxima_simple(
            &maxima_coords_filtered,
            &maxima_vals,
            min_distance,
            cell_ref.as_ref(),
            (nx, ny, nz),
            false, // non-periodic
        )
    } else {
        maxima_coords_filtered.clone() // No merging  
    };
    
    // Convert to NumPy arrays
    let n_raw = maxima_coords_filtered.len();
    let raw_arr = Array2::<usize>::from_shape_vec((n_raw, 3), maxima_coords_filtered.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_pbc_arr = Array2::<usize>::from_shape_vec((merged_pbc_coords.len(), 3), merged_pbc_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    let merged_cart_arr = Array2::<usize>::from_shape_vec((merged_cart_coords.len(), 3), merged_cart_coords.iter().flat_map(|c| c.iter()).cloned().collect()).unwrap();
    
    Ok((
        PyArray2::from_owned_array(py, merged_pbc_arr).to_owned(),
        PyArray2::from_owned_array(py, merged_cart_arr).to_owned(),
        PyArray2::from_owned_array(py, raw_arr).to_owned(),
    ))
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
    m.add_function(wrap_pyfunction!(detect_maxima_h_morph, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_h_morph_opening, m)?)?;
    m.add_function(wrap_pyfunction!(move_maxima_to_peaks, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_h_morph2, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_h_morph_simple, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_h_morph_no_spacing, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_adaptive_geometry, m)?)?;
    m.add_function(wrap_pyfunction!(detect_maxima_selective_smoothing, m)?)?;
    m.add("__all__", ["dgrid_from_atoms", "egrid_from_atoms", "window_connections", "periodic_gaussian_filter", "seeded_watershed_bucket", "seeded_watershed", "seeded_watershed_masked", "detect_maxima_h_morph", "detect_maxima_h_morph_opening", "move_maxima_to_peaks", "detect_maxima_h_morph2", "detect_maxima_h_morph_simple", "detect_maxima_h_morph_no_spacing", "detect_maxima_adaptive_geometry", "detect_maxima_selective_smoothing"])?;

    // By returning Ok(()), the module initialization is successful.
    Ok(())
} 