use std::io::{BufRead, BufReader};
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::time::Instant;
use std::collections::{HashMap, HashSet, BTreeSet};
use std::collections::hash_map::Entry;
use std::collections::hash_map::DefaultHasher;
use std::num::Wrapping;
use std::hash::{Hash, Hasher};
// use std::cmp::{max};
use rand::Rng;

/// Return structure holding final stats
#[derive(Debug)]
pub struct SignificanceSummary {
    pub num_significant: usize,
    pub num_not_significant: usize,

    pub mean_cov_significant: f64,
    pub std_cov_significant: f64,

    pub mean_cov_not_significant: f64,
    pub std_cov_not_significant: f64,
}

#[derive(Debug, Clone)]
pub struct Pair<T, U> {
    pub value0: T,
    pub value1: U,
}

impl<T, U> Pair<T, U> {
    pub fn new(value0: T, value1: U) -> Self {
        Pair { value0, value1 }
    }
}

#[derive(Debug, Clone)]
pub struct Triple<T, U, V> {
    pub value0: T,
    pub value1: U,
    pub value2: V,
}

impl<T, U, V> Triple<T, U, V> {
    pub fn new(value0: T, value1: U, value2: V) -> Self {
        Triple { value0, value1, value2 }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Quartet<T, U, V, W> {
    pub value0: T,
    pub value1: U,
    pub value2: V,
    pub value3: W,
}

impl<T, U, V, W> Quartet<T, U, V, W> {
    pub fn new(value0: T, value1: U, value2: V, value3: W) -> Self {
        Quartet { value0, value1, value2, value3 }
    }
}


pub fn get_lines(file_name: &str, size: usize) -> Vec<Vec<String>> {
    let mut records = Vec::new();
    if let Ok(file) = File::open(file_name) {
        let reader = BufReader::new(file);
        for line_result in reader.lines() {
            if let Ok(line) = line_result {
                let fields: Vec<String> = line
                    .split(',')
                    .map(|s| s.to_string())
                    .collect();
                
                // check for missing values '?'
                let mut has_missing = false;
                for s in &fields {
                    if s.contains('?') {
                        has_missing = true;
                        break;
                    }
                }
                if has_missing {
                    continue;
                }

                // The Java code changes 'e' -> "1", 'p' -> "0" if they appear at fields[size].
                // "fields[SIZE]" is the (SIZE)-th index (0-based) in the array.
                // Here we do the same logic:
                let mut fields_copy = fields.clone();
                if fields_copy[size] == "e" {
                    fields_copy[size] = "1".to_string();
                } else if fields_copy[size] == "p" {
                    fields_copy[size] = "0".to_string();
                } else {
                    // smth = 0 in Java. We ignore it or handle it as needed.
                    // For now, we do nothing. 
                }

                records.push(fields_copy);
            }
        }
    }
    records
}


pub fn create_index(
    records: &Vec<Vec<String>>,
    size: usize,
    num_labels: usize,
    index_holder: &mut HashMap<String, i32>,
    lengths: &mut Vec<i32>,
    enum_recs: &mut Vec<Vec<i32>>,
    binary_labels: &mut Vec<Vec<i32>>,
) {
    let mut first_line = true;
    for record in records {
        if first_line {
            first_line = false;
            // In the Java code, we skip the first record if it might be a header
            continue;
        }
        let mut enum_rec = Vec::with_capacity(size);
        let mut bin_label = vec![0; num_labels];
        
        for (i, val) in record.iter().enumerate() {
            if i >= size {
                // put into binary_labels
                bin_label[i - size] = val.parse::<i32>().unwrap_or(0);
            } else {
                // build key
                let key = format!("{} {}", i, val);
                if !index_holder.contains_key(&key) {
                    lengths[i] += 1;
                    index_holder.insert(key.clone(), lengths[i]);
                }
                let mapped_value = *index_holder.get(&key).unwrap();
                enum_rec.push(mapped_value);
            }
        }
        enum_recs.push(enum_rec);
        binary_labels.push(bin_label);
    }
}


pub fn get_mover(lengths: &Vec<i32>, size: usize) -> Vec<i32> {
    let mut bit_size = vec![0; size];
    for i in 1..size {
        // bit_size[i] = bit_size[i-1] + (int)(Math.log(lengths[i-1])/Math.log(2)) + 1;
        let prev = bit_size[i - 1];
        let log2_val = (lengths[i - 1] as f64).log2();
        bit_size[i] = prev + log2_val as i32 + 1;
    }
    bit_size
}


pub fn next_mask(current: i64) -> i64 {
    let c = current & -current;
    let r = current + c;
    (((r ^ current) >> 2) / c) | r
}


// pub fn aggregate(
//     mov: &Vec<i32>,
//     _lengths: &Vec<i32>,
//     _index_holder: &HashMap<String, i32>, // not deeply used in your snippet
//     map: &mut HashMap<i64, Vec<i64>>,
//     enum_recs: &Vec<Vec<i32>>,
//     binary_labels: &Vec<Vec<i32>>,
//     covers: &mut HashMap<i64, i64>,
//     size: usize,
//     num_labels: usize,
//     threshold: f64,
// ) {
//     let mut agg: HashMap<Vec<i32>, Vec<i64>> = HashMap::new();
//     let hash: i64 = 1;
//     let mut count = 0i64;
//     let large_prime: i64 = 2305843009213693951; // 2^61 - 1

//     let mut bases: HashMap<Vec<i32>, i64> = HashMap::new();

//     for list_key in enum_recs {
//         if !agg.contains_key(list_key) {
//             // create new label array
//             let mut labels = vec![0i64; num_labels];
//             for y in 0..num_labels {
//                 // in Java: (1L << 32) + binaryLabels[count][y]
//                 labels[y] = (1i64 << 32) + (binary_labels[count as usize][y] as i64);
//             }
//             agg.insert(list_key.clone(), labels);

//             // compute base
//             let mut base: i64 = 0;
//             for (j, val) in list_key.iter().enumerate().take(size) {
//                 base += (*val as i64) << mov[j];
//             }
//             bases.insert(list_key.clone(), base);
//         } else {
//             // update existing label array
//             let existing = agg.get_mut(list_key).unwrap();
//             for y in 0..num_labels {
//                 existing[y] += (1i64 << 32) + (binary_labels[count as usize][y] as i64);
//             }
//         }
//         count += 1;
//     }

//     // println!("Agg size: {}", agg.len());

//     for i in 0..=size {
//         let mut cov: HashMap<i64, i64> = HashMap::new();
//         let mut stat: HashMap<i64, Vec<i64>> = HashMap::new();

//         for list_key in agg.keys() {
//             let base = bases[list_key];
//             let labels = &agg[list_key];

//             if i == 0 {
//                 // everything merges into key 0
//                 let l = 0i64;
//                 map.entry(l).or_insert(vec![0; num_labels]);
//                 covers.entry(l).or_insert(hash);
                
//                 let entry_map = map.get_mut(&l).unwrap();
//                 for y in 0..num_labels {
//                     entry_map[y] += labels[y];
//                 }

//                 let cv = covers.get_mut(&l).unwrap();
//                 *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64; 

//             } else if i == size {
//                 // only add if #inc >= threshold
//                 // from Java: (stat.get(l)[0] >> 32) is the total
//                 // (stat.get(l)[0] % (1 << 32)) is inc but the code is a bit different
//                 // Actually: getRate(long) => inc / tot
//                 // We check if tot >= threshold or inc >= threshold?
//                 // The Java snippet does: if (double)(agg.get(listKey)[0] >> 32) >= threshold
//                 // That means if total >= threshold. 
//                 let tot = (labels[0] >> 32) as f64;
//                 if tot >= threshold {
//                     map.insert(base, labels.clone());
//                     let cval = (hash * (base + 31)) % large_prime;
//                     covers.insert(base, cval);
//                 }

//             } else {
//                 // loop over subsets of size i
//                 // mask from (1 << i) - 1 up to ...
//                 let start_mask: i64 = ((1i64 << i) - 1) << (size - i);
//                 // The Java code: for (long mask = (1L << i) - 1; mask <= ((1L << i) - 1) << (SIZE - i); ...)
//                 // We must be careful about how we iterate in Rust because we do next_mask.
//                 // So we do something like:

//                 let mut mask = (1i64 << i) - 1;
//                 while mask <= start_mask {
//                     let mut l = 0i64;
//                     for j in 0..size {
//                         let ind = 1i64 << j;
//                         if (ind & mask) == ind {
//                             l += (list_key[j] as i64) << mov[j];
//                         }
//                     }

//                     // now update stat and cov
//                     match stat.entry(l) {
//                         Entry::Vacant(e) => {
//                             e.insert(labels.clone());
//                             cov.insert(l, hash);
//                             // update covers
//                             let cv = cov.get_mut(&l).unwrap();
//                             // *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64;
//                             *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

//                             // check threshold
//                             let tot = (labels[0] >> 32) as f64;
//                             if tot >= threshold {
//                                 map.insert(l, labels.clone());
//                                 covers.insert(l, *cv);
//                             }
//                         }
//                         Entry::Occupied(mut e) => {
//                             // If already in map but not in the global map
//                             // or partially in the global map
//                             let e_val = e.get_mut();
//                             for y in 0..num_labels {
//                                 e_val[y] += labels[y];
//                             }

//                             let cv = cov.get_mut(&l).unwrap();
//                             // *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64;
//                             *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

//                             let tot = (e_val[0] >> 32) as f64;
//                             if tot >= threshold {
//                                 // update global
//                                 map.insert(l, e_val.clone());
//                                 covers.insert(l, *cv);
//                             }
//                         }
//                     }

//                     // proceed to next mask
//                     mask = next_mask(mask);
//                 }
//             }
//         }
//     }
// }


pub fn aggregate(
    mov: &Vec<i32>,
    _lengths: &Vec<i32>,            // not deeply used in your snippet
    _index_holder: &HashMap<String, i32>,
    map: &mut HashMap<i64, Vec<i64>>,    // final global map to be populated
    enum_recs: &Vec<Vec<i32>>,
    binary_labels: &Vec<Vec<i32>>,
    covers: &mut HashMap<i64, i64>, // final global covers to be populated
    size: usize,
    num_labels: usize,
    threshold: f64,
) {
    let mut agg: HashMap<Vec<i32>, Vec<i64>> = HashMap::new();
    let mut bases: HashMap<Vec<i32>, i64> = HashMap::new();

    let hash: i64 = 1;
    let large_prime: i64 = 2305843009213693951; // 2^61 - 1

    // 1) BUILD agg and bases (sequential)
    let mut count = 0i64;
    for list_key in enum_recs {
        if !agg.contains_key(list_key) {
            let mut labels = vec![0i64; num_labels];
            for y in 0..num_labels {
                // (1 << 32) + ...
                labels[y] = (1i64 << 32) + (binary_labels[count as usize][y] as i64);
            }
            agg.insert(list_key.clone(), labels);

            // compute base
            let mut base: i64 = 0;
            for (j, val) in list_key.iter().enumerate().take(size) {
                base += (*val as i64) << mov[j];
            }
            bases.insert(list_key.clone(), base);
        } else {
            // update existing
            let existing = agg.get_mut(list_key).unwrap();
            for y in 0..num_labels {
                existing[y] += (1i64 << 32) + (binary_labels[count as usize][y] as i64);
            }
        }
        count += 1;
    }

    println!("Agg size: {}", agg.len());

    // 2) PARALLELIZE the loop over i in 0..=size
    //    For each i, we build a local stat/cov, do the logic, and return partial results.

    // We'll define a closure that processes one "i" value and returns partial maps
    //    partial_map: HashMap<i64, Vec<i64>>
    //    partial_covers: HashMap<i64, i64>

    let partial_results: Vec<(HashMap<i64, Vec<i64>>, HashMap<i64, i64>)> =
        (0..=size).into_par_iter()
        .map(|i| {
            // Each parallel task gets its own local map
            let mut local_map: HashMap<i64, Vec<i64>> = HashMap::new();
            let mut local_covers: HashMap<i64, i64> = HashMap::new();

            let mut stat: HashMap<i64, Vec<i64>> = HashMap::new();
            let mut cov: HashMap<i64, i64> = HashMap::new();

            // "agg" and "bases" are shared read-only
            for (list_key, labels) in &agg {
                let base = *bases.get(list_key).unwrap();

                if i == 0 {
                    // merges everything into key=0
                    let l = 0i64;
                    local_map.entry(l).or_insert(vec![0; num_labels]);

                    local_covers.entry(l).or_insert(hash);

                    // update local_map
                    let entry_map = local_map.get_mut(&l).unwrap();
                    for y in 0..num_labels {
                        entry_map[y] += labels[y];
                    }
                    // update local_covers
                    let cv = local_covers.get_mut(&l).unwrap();
                    // e.g. wrapping multiply
                    *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

                } else if i == size {
                    // only add if total >= threshold
                    let tot = (labels[0] >> 32) as f64;
                    if tot >= threshold {
                        local_map.insert(base, labels.clone());
                        let cval = (hash * (base + 31)) % large_prime;
                        local_covers.insert(base, cval);
                    }
                } else {
                    // subsets of size i
                    let start_mask: i64 = ((1i64 << i) - 1) << (size - i);
                    let mut mask = (1i64 << i) - 1;
                    while mask <= start_mask {
                        let mut l = 0i64;
                        for j in 0..size {
                            let ind = 1i64 << j;
                            if (ind & mask) == ind {
                                l += (list_key[j] as i64) << mov[j];
                            }
                        }

                        match stat.entry(l) {
                            Entry::Vacant(e) => {
                                e.insert(labels.clone());
                                cov.insert(l, hash);

                                let cv = cov.get_mut(&l).unwrap();
                                *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

                                let tot = (labels[0] >> 32) as f64;
                                if tot >= threshold {
                                    local_map.insert(l, labels.clone());
                                    local_covers.insert(l, *cv);
                                }
                            }
                            Entry::Occupied(mut e) => {
                                let e_val = e.get_mut();
                                for y in 0..num_labels {
                                    e_val[y] += labels[y];
                                }

                                let cv = cov.get_mut(&l).unwrap();
                                *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

                                let tot = (e_val[0] >> 32) as f64;
                                if tot >= threshold {
                                    local_map.insert(l, e_val.clone());
                                    local_covers.insert(l, *cv);
                                }
                            }
                        }

                        mask = next_mask(mask);
                    }
                }
            } // end for (list_key, labels) in &agg

            // Return partial results for this i
            (local_map, local_covers)
        })
        .collect();

    // 3) MERGE partial results into the "global" map/covers

    for (pm, pc) in partial_results {
        // Merge pm into map
        for (k, v) in pm {
            // If the key does not exist in global map, insert it
            match map.entry(k) {
                Entry::Vacant(e) => {
                    e.insert(v);
                }
                Entry::Occupied(mut oe) => {
                    // If it exists, we can merge by summation or do something else
                    let existing_vec = oe.get_mut();
                    for (idx, val) in v.iter().enumerate() {
                        existing_vec[idx] += *val;
                    }
                }
            }
        }

        // Merge pc into covers
        for (k, v) in pc {
            match covers.entry(k) {
                Entry::Vacant(e) => {
                    e.insert(v);
                }
                Entry::Occupied(mut oe) => {
                    // if needed, combine coverage 
                    // e.g. keep a product or something
                    // Below we just overwrite, but you can define your logic
                    let existing_val = oe.get_mut();
                    *existing_val = (*existing_val).wrapping_add(v) % large_prime;
                }
            }
        }
    }

    // done
    // println!("Done parallel loop for i in 0..=size");
}


pub fn get_rate(l2: i64) -> f64 {
    let inc = (l2 & 0xFFFF_FFFF) as f64; // l2 % (1 << 32)
    let tot = (l2 >> 32) as f64;
    if tot == 0.0 {
        0.0
    } else {
        inc / tot 
    }
}


/// Groups covers so that the covering (i64) is the key
/// and the value is a `Vec<i64>` of populations having that cover.
pub fn group_covers(
    covers: &HashMap<i64, i64>
) -> HashMap<i64, Vec<i64>> {
    let mut gc: HashMap<i64, Vec<i64>> = HashMap::new();
    // covers: key = population, value = cover
    for (pop, cover_val) in covers {
        let new_key = *cover_val;
        gc.entry(new_key)
          .and_modify(|pop_vec| {
              pop_vec.push(*pop);
          })
          .or_insert({
              let mut v = Vec::new();
              v.push(*pop);
              v
          });
    }
    gc
}


/// Retrieve the bitwise indexes for a population `key`
pub fn get_indexes(key: i64, mover: &Vec<i32>, size: usize) -> Vec<i64> {
    let mut ret = vec![0i64; size];
    for j in 0..(size - 1) {
        let width = mover[j + 1] - mover[j];
        ret[j] = (key >> mover[j]) & ((1 << width) - 1);
    }
    ret[size - 1] = key >> mover[size - 1];
    ret
}


pub fn check_sp(
    c1: i64,
    c2: i64,
    j: usize,
    y: usize,
    stats: &HashMap<i64, Vec<i64>>,
    mov: &Vec<i32>,
    lengths: &Vec<i32>,
    covers: &HashMap<i64, i64>,
) -> Pair<bool, i64> {
    let large_prime: i64 = 2305843009213693951; // 2^61 - 1
    let mut hash_acc = 1i64;

    // let mut separators_involved = Vec::new();

    // if stats does not contain c1 or c2, return false
    if !stats.contains_key(&c1) || !stats.contains_key(&c2) {
        return Pair::new(false, hash_acc);
    }
    let cover_c1 = *covers.get(&c1).unwrap_or(&0);
    hash_acc = hash_acc.wrapping_mul(cover_c1.wrapping_add(31)) % large_prime;

    let cover_c2 = *covers.get(&c2).unwrap_or(&0);
    hash_acc = hash_acc.wrapping_mul(cover_c2.wrapping_add(31)) % large_prime;

    let c1r = get_rate(stats[&c1][y]);
    let bits_c1r = c1r.to_bits() as i64;
    hash_acc = hash_acc.wrapping_mul(bits_c1r.wrapping_add(31)) % large_prime;

    let c2r = get_rate(stats[&c2][y]);
    let bits_c2r = c2r.to_bits() as i64;
    hash_acc = hash_acc.wrapping_mul(bits_c2r.wrapping_add(31)) % large_prime;

    if c1r > c2r {
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return Pair::new(false, hash_acc);
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr > c2xr {
                return Pair::new(false, hash_acc);
            }
            let cover_c1x = *covers.get(&c1x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c1x.wrapping_add(31)) % large_prime;

            let cover_c2x = *covers.get(&c2x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c2x.wrapping_add(31)) % large_prime;

            let bits_c1xr = c1xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c1xr.wrapping_add(31)) % large_prime;

            let bits_c2xr = c2xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c2xr.wrapping_add(31)) % large_prime;
        }
    } else if c1r < c2r {
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return Pair::new(false, hash_acc);
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr < c2xr {
                return Pair::new(false, hash_acc);
            }

            let cover_c1x = *covers.get(&c1x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c1x.wrapping_add(31)) % large_prime;

            let cover_c2x = *covers.get(&c2x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c2x.wrapping_add(31)) % large_prime;

            let bits_c1xr = c1xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c1xr.wrapping_add(31)) % large_prime;

            let bits_c2xr = c2xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c2xr.wrapping_add(31)) % large_prime;
        }
    } else {
        // c1r == c2r
        let mut c1gc2 = 0;
        let mut c2gc1 = 0;
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return Pair::new(false, hash_acc);
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr < c2xr {
                c2gc1 += 1;
            } else if c1xr > c2xr {
                c1gc2 += 1;
            } else {
                // c1xr == c2xr
                return Pair::new(false, hash_acc);
            }

            let cover_c1x = *covers.get(&c1x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c1x.wrapping_add(31)) % large_prime;

            let cover_c2x = *covers.get(&c2x).unwrap_or(&0);
            hash_acc = hash_acc.wrapping_mul(cover_c2x.wrapping_add(31)) % large_prime;

            let bits_c1xr = c1xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c1xr.wrapping_add(31)) % large_prime;

            let bits_c2xr = c2xr.to_bits() as i64;
            hash_acc = hash_acc.wrapping_mul(bits_c2xr.wrapping_add(31)) % large_prime;

        }
        // (c1gc2 == 0 || c2gc1 == 0)
        return Pair::new(c1gc2 == 0 || c2gc1 == 0, hash_acc);
    }

    Pair::new(true, hash_acc)
}


pub fn check_sp_sign(
    c1: i64,
    c2: i64,
    j: usize,
    y: usize,
    stats: &HashMap<i64, Vec<i64>>,
    mov: &Vec<i32>,
    lengths: &Vec<i32>
) -> bool {

    // if stats does not contain c1 or c2, return false
    if !stats.contains_key(&c1) || !stats.contains_key(&c2) {
        return false;
    }
    let c1r = get_rate(stats[&c1][y]);

    let c2r = get_rate(stats[&c2][y]);

    if c1r > c2r {
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return false;
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr > c2xr {
                return false;
            }
        }
    } else if c1r < c2r {
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return false;
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr < c2xr {
                return false;
            }
        }
    } else {
        // c1r == c2r
        let mut c1gc2 = 0;
        let mut c2gc1 = 0;
        for k in 1..=lengths[j] {
            let c1x = c1 + ((k as i64) << mov[j]);
            let c2x = c2 + ((k as i64) << mov[j]);
            if !stats.contains_key(&c1x) || !stats.contains_key(&c2x) {
                return false;
            }
            let c1xr = get_rate(stats[&c1x][y]);
            let c2xr = get_rate(stats[&c2x][y]);

            if c1xr < c2xr {
                c2gc1 += 1;
            } else if c1xr > c2xr {
                c1gc2 += 1;
            } else {
                // c1xr == c2xr
                return false;
            }

        }
        // (c1gc2 == 0 || c2gc1 == 0)
        return c1gc2 == 0 || c2gc1 == 0;
    }

    true
}


/// Returns a Pair of:
/// - A vector of indices that are zero for the best candidate
/// - The population key with the fewest zeros
pub fn get_least_stars(
    list: &Vec<i64>,
    mover: &Vec<i32>,
    size: usize,
) -> Pair<Vec<usize>, i64> {
    let mut least = size + 1;
    let mut key = 0i64;
    let mut zeros = Vec::new();

    // Iterate over all population keys in `list`
    for &l in list.iter() {
        let mut num_zeros = 0;
        let mut nzeros = Vec::new();

        // For each attribute except the last one, we measure the bits
        for j in 0..(size - 1) {
            let width = mover[j + 1] - mover[j];
            let val = (l >> mover[j]) & ((1 << width) - 1);
            if val == 0 {
                num_zeros += 1;
                nzeros.push(j);
            }
        }

        // Check the last chunk
        let val_last = l >> mover[size - 1];
        if val_last == 0 {
            num_zeros += 1;
            nzeros.push(size - 1);
        }

        if num_zeros < least {
            least = num_zeros;
            key = l;
            zeros = nzeros;
        }
    }

    Pair::new(zeros, key)
}


// pub fn find_sp(
//     aggregations: &HashMap<i64, Vec<i64>>,
//     grouped_covers: &HashMap<i64, Vec<i64>>,   // changed here
//     lengths: &Vec<i32>,
//     mover: &Vec<i32>,
//     covers: &HashMap<i64, i64>,
//     size: usize,
//     num_labels: usize,
// ) -> Pair<
//     HashSet<Quartet<i64, i64, usize, usize>>,
//     HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>>
// >
// {
//     let mut all_sps: HashSet<Quartet<i64, i64, usize, usize>> = HashSet::new();
//     let mut infos: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();
//     let mut infos_aux: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();

//     // groupedCovers: key = cover value, value = Vec<i64> of populations
//     for (_cover_val, pop_list) in grouped_covers {
//         // getLeastStars now uses &Vec<i64>
//         let stars = get_least_stars(pop_list, mover, size);
//         let l1 = stars.value1;
//         let zeros = stars.value0;

//         if zeros.len() < 2 {
//             continue;
//         }

//         // for each attribute i in zeros
//         for &i in zeros.iter() {
//             // i1 in [1..lengths[i]]
//             for i1 in 1..lengths[i] {
//                 let c1 = l1 + ((i1 as i64) << mover[i]);
//                 if !aggregations.contains_key(&c1) {
//                     continue;
//                 }

//                 for i2 in (i1 + 1)..=lengths[i] {
//                     let c2 = l1 + ((i2 as i64) << mover[i]);
//                     if !aggregations.contains_key(&c2) {
//                         continue;
//                     }

//                     // for each potential separator j
//                     for &j in zeros.iter() {
//                         if j == i {
//                             continue;
//                         }

//                         for y in 0..num_labels {
//                             let p = Quartet::new(c1, c2, j, y);
//                             if all_sps.contains(&p) {
//                                 continue;
//                             }
//                             let res = check_sp(c1, c2, j, y, aggregations, mover, lengths, covers);
//                             let is_sp = res.value0;
//                             let info = res.value1;

//                             if is_sp && !infos_aux.contains_key(&info) {
//                                 // create a new "info bucket"
//                                 all_sps.insert(p.clone());
//                                 let mut new_set = HashSet::new();
//                                 new_set.insert(p.clone());
//                                 infos_aux.insert(info.clone(), new_set.clone());
//                                 infos.insert(info.clone(), new_set);

//                                 // find siblings
//                                 let c1cov = covers.get(&c1).unwrap();
//                                 let c2cov = covers.get(&c2).unwrap();

//                                 let c1cov_list = grouped_covers.get(c1cov).unwrap();
//                                 let c2cov_list = grouped_covers.get(c2cov).unwrap();

//                                 // examine all pairs in c1cov_list x c2cov_list
//                                 for &t1 in c1cov_list.iter() {
//                                     let t1_ind = get_indexes(t1, mover, size);

//                                     for &t2 in c2cov_list.iter() {
//                                         let t2_ind = get_indexes(t2, mover, size);
//                                         let mut diff = 0;
//                                         for idx in 0..size {
//                                             if t1_ind[idx] != t2_ind[idx] {
//                                                 diff += 1;
//                                             }
//                                         }
//                                         if diff == 1 {
//                                             let small = if t1 < t2 { t1 } else { t2 };
//                                             let big = if t1 < t2 { t2 } else { t1 };

//                                             let p_prime = Quartet::new(small, big, j, y);
//                                             all_sps.insert(p_prime.clone());

//                                             let info_set = infos.get_mut(&info).unwrap();
//                                             info_set.insert(p_prime.clone());
//                                             let info_aux_set = infos_aux.get_mut(&info).unwrap();
//                                             info_aux_set.insert(p_prime);
//                                         }
//                                     }
//                                 }
//                             } else if is_sp {
//                                 // If is_sp and info already present, copy in the new p or merge
//                                 if let Some(old_set) = infos_aux.get(&info) {
//                                     for p_prime in old_set.iter() {
//                                         let p_hat = Quartet::new(p_prime.value0, p_prime.value1, j, y);
//                                         all_sps.insert(p_hat.clone());
//                                         if let Some(infos_set) = infos.get_mut(&info) {
//                                             infos_set.insert(p_hat);
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     Pair::new(all_sps, infos)
// }


// -------------------------------------------------------------------
// Parallel version of find_sp
// -------------------------------------------------------------------
pub fn find_sp(
    aggregations: &HashMap<i64, Vec<i64>>,
    grouped_covers: &HashMap<i64, Vec<i64>>,
    lengths: &Vec<i32>,
    mover: &Vec<i32>,
    covers: &HashMap<i64, i64>,
    size: usize,
    num_labels: usize,
) -> Pair<
    HashSet<Quartet<i64, i64, usize, usize>>,
    HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>>
>
{
    // We will build partial results for each cover_val, then merge them.
    // Each partial result is (HashSet<Quartet<...>>, HashMap<i64, HashSet<Quartet<...>>>).

    let partial_results: Vec<(
        HashSet<Quartet<i64, i64, usize, usize>>,
        HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>>
    )> = grouped_covers
        .par_iter() // parallel iteration over (cover_val, pop_list)
        .map(|(_cover_val, pop_list)| {
            // Build local sets for this single cover_val
            let mut local_sps: HashSet<Quartet<i64, i64, usize, usize>> = HashSet::new();
            let mut local_infos: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();

            // We also need a local "aux" if your logic uses it:
            let mut local_infos_aux: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();

            // replicate the logic from your original loop body:
            let stars = get_least_stars(pop_list, mover, size);
            let l1 = stars.value1;
            let zeros = stars.value0;

            if zeros.len() < 2 {
                // Return empty partial sets
                return (local_sps, local_infos);
            }

            for &i in zeros.iter() {
                for i1 in 1..lengths[i] {
                    let c1 = l1 + ((i1 as i64) << mover[i]);
                    if !aggregations.contains_key(&c1) {
                        continue;
                    }

                    for i2 in (i1 + 1)..=lengths[i] {
                        let c2 = l1 + ((i2 as i64) << mover[i]);
                        if !aggregations.contains_key(&c2) {
                            continue;
                        }

                        // for each potential separator j
                        for &j in zeros.iter() {
                            if j == i {
                                continue;
                            }
                            for y in 0..num_labels {
                                let p = Quartet::new(c1, c2, j, y);
                                if local_sps.contains(&p) {
                                    continue;
                                }
                                let res = check_sp(c1, c2, j, y, aggregations, mover, lengths, covers);
                                let is_sp = res.value0;
                                let info_hash = res.value1; // i64

                                if is_sp && !local_infos_aux.contains_key(&info_hash) {
                                    // create a new "info bucket"
                                    local_sps.insert(p.clone());
                                    let mut new_set = HashSet::new();
                                    new_set.insert(p.clone());
                                    local_infos_aux.insert(info_hash, new_set.clone());
                                    local_infos.insert(info_hash, new_set);

                                    // find siblings
                                    let c1cov = covers.get(&c1).unwrap();
                                    let c2cov = covers.get(&c2).unwrap();

                                    let c1cov_list = grouped_covers.get(c1cov).unwrap();
                                    let c2cov_list = grouped_covers.get(c2cov).unwrap();

                                    for &t1 in c1cov_list.iter() {
                                        let t1_ind = get_indexes(t1, mover, size);
                                        for &t2 in c2cov_list.iter() {
                                            let t2_ind = get_indexes(t2, mover, size);
                                            let mut diff = 0;
                                            for idx in 0..size {
                                                if t1_ind[idx] != t2_ind[idx] {
                                                    diff += 1;
                                                }
                                            }
                                            if diff == 1 {
                                                let small = if t1 < t2 { t1 } else { t2 };
                                                let big   = if t1 < t2 { t2 } else { t1 };

                                                let p_prime = Quartet::new(small, big, j, y);
                                                local_sps.insert(p_prime.clone());

                                                // update local_infos & local_infos_aux
                                                let info_set = local_infos.get_mut(&info_hash).unwrap();
                                                info_set.insert(p_prime.clone());
                                                let info_aux_set = local_infos_aux.get_mut(&info_hash).unwrap();
                                                info_aux_set.insert(p_prime);
                                            }
                                        }
                                    }
                                } else if is_sp {
                                    // If is_sp and info already present, copy in the new p or merge
                                    if let Some(old_set) = local_infos_aux.get(&info_hash) {
                                        for p_prime in old_set.iter() {
                                            let p_hat = Quartet::new(p_prime.value0, p_prime.value1, j, y);
                                            local_sps.insert(p_hat.clone());
                                            if let Some(infos_set) = local_infos.get_mut(&info_hash) {
                                                infos_set.insert(p_hat);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Return partial results for this single cover_val
            (local_sps, local_infos)
        })
        .collect(); // gather all partial results into a Vec

    // Now, merge partial results
    let mut all_sps: HashSet<Quartet<i64, i64, usize, usize>> = HashSet::new();
    let mut infos: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();

    for (pset, imap) in partial_results {
        // Merge pset into all_sps
        all_sps.extend(pset);

        // Merge imap into infos
        for (info_hash, qset) in imap {
            match infos.entry(info_hash) {
                Entry::Vacant(e) => {
                    e.insert(qset);
                }
                Entry::Occupied(mut oe) => {
                    let existing = oe.get_mut();
                    existing.extend(qset);
                }
            }
        }
    }

    // Return final merged result
    Pair::new(all_sps, infos)
}

/// Merged function that both decodes the 64-bit population into attribute codes
/// and counts which columns are zero-coded (stars).
/// Returns (star_columns, star_count).
///
/// - `pop`: the 64-bit population key
/// - `size`: number of attributes
/// - `mover`: bit offsets for each attribute
fn star_info_for_pop(pop: i64, size: usize, mover: &Vec<i32>) -> (HashSet<usize>, i32) {
    let mut star_cols = HashSet::new();
    let mut star_count = 0;

    // Decode attribute codes in one pass, check if zero => star
    for j in 0..(size - 1) {
        let width = mover[j + 1] - mover[j];
        let code_val = (pop >> mover[j]) & ((1 << width) - 1);
        if code_val == 0 {
            star_cols.insert(j);
            star_count += 1;
        }
    }
    // last attribute:
    let code_val = pop >> mover[size - 1];
    if code_val == 0 {
        star_cols.insert(size - 1);
        star_count += 1;
    }

    (star_cols, star_count)
}

/// Finds the shared star attributes among a single group of redundant SP:
/// a HashSet of `Quartet<i64, i64, usize, usize>`. We only look at `value0` in each paradox.
/// Also finds the paradox with the fewest star columns and the one with the most star columns.
///
/// Returns a triple:
/// (shared_star_columns, min_star_paradox, max_star_paradox).
pub fn find_shared_star_attribute_sibling(
    group: &HashSet<Quartet<i64, i64, usize, usize>>,
    size: usize,
    mover: &Vec<i32>,
) -> (
    Vec<usize>,                                 // shared star columns
    Option<Quartet<i64, i64, usize, usize>>,    // paradox with fewest stars
    Option<Quartet<i64, i64, usize, usize>>     // paradox with most stars
) {
    if group.is_empty() {
        return (vec![], None, None);
    }

    // Convert set to iterator
    let mut iter = group.iter();

    // Take first paradox as baseline
    let first_sp = iter.next().unwrap();
    let (mut shared_cols, mut min_count) = star_info_for_pop(first_sp.value0, size, mover);
    let mut max_count = min_count;
    let mut min_sp: Option<Quartet<i64, i64, usize, usize>> = Some(first_sp.clone());
    let mut max_sp: Option<Quartet<i64, i64, usize, usize>> = Some(first_sp.clone());

    let mut division_cols = HashSet::new();
    division_cols.insert(first_sp.value2);

    // For each subsequent paradox, intersect star cols and update min/max
    for sp in iter {
        let (sp_cols, sp_count) = star_info_for_pop(sp.value0, size, mover);

        // Union the star columns
        shared_cols.extend(&sp_cols);

        // Update min star
        if sp_count < min_count {
            min_count = sp_count;
            min_sp = Some(sp.clone());
        }
        // Update max star
        if sp_count > max_count {
            max_count = sp_count;
            max_sp = Some(sp.clone());
        }

        division_cols.insert(sp.value2);

        if shared_cols.is_empty() {
            // No need to continue
            break;
        }
    }

    // Convert from HashSet to Vec
    shared_cols.retain(|col| !division_cols.contains(col));
    let sibling_vec: Vec<usize> = shared_cols.into_iter().collect();
    (sibling_vec, min_sp, max_sp)
}


pub fn find_shared_star_attribute_division(
    group: &HashSet<Quartet<i64, i64, usize, usize>>,
    size: usize,
    mover: &Vec<i32>,
) -> (
    Vec<usize>,                                 // division_vec (the union of value2 from all paradoxes)
    Option<Quartet<i64, i64, usize, usize>>,    // sp_diff1
    Option<Quartet<i64, i64, usize, usize>>     // sp_diff2
) {
    if group.is_empty() {
        return (vec![], None, None);
    }

    // Convert set to iterator
    let mut iter = group.iter();

    // Take first paradox as baseline
    let first_sp = iter.next().unwrap();
    let (mut shared_cols, mut min_count) = star_info_for_pop(first_sp.value0, size, mover);
    let mut min_sp: Option<Quartet<i64, i64, usize, usize>> = Some(first_sp.clone());
    let mut max_sp: Option<Quartet<i64, i64, usize, usize>> = Some(first_sp.clone());

    let mut division_cols = HashSet::new();
    division_cols.insert(first_sp.value2);

    // for each subsequent paradox, union star cols and union the value2
    for sp in iter {
        let (sp_cols, sp_count) = star_info_for_pop(sp.value0, size, mover);

        // Union the star columns
        shared_cols.extend(&sp_cols);

        division_cols.insert(sp.value2);

        // If we have not picked sp_diff2 yet, we check if sp.value2 differs from first_sp.value2
        if sp.value2 != min_sp.clone().unwrap().value2 {
            max_sp = Some(sp.clone());
        }
    }

    // division_vec => all distinct value2 codes
    let shared_vec: Vec<usize> = division_cols.into_iter().collect();
    // let shared_vec: Vec<usize> = shared_cols.into_iter().collect();
    (shared_vec, min_sp, max_sp)
}


/// If record_codes matches `pop` under star logic, returns true.
/// Star logic means any code=0 in `pop` is do-not-care for that attribute.
fn record_matches_pop_with_stars(
    record_codes: &Vec<i32>,
    pop: i64,
    size: usize,
    mover: &Vec<i32>,
) -> bool {
    for j in 0..(size - 1) {
        let width = mover[j + 1] - mover[j];
        let code_val = (pop >> mover[j]) & ((1 << width) - 1);
        if code_val != 0 {
            // must match exactly
            if code_val as i32 != record_codes[j] {
                return false;
            }
        }
    }
    // last attribute
    let code_val_last = pop >> mover[size - 1];
    if code_val_last != 0 {
        if code_val_last as i32 != record_codes[size - 1] {
            return false;
        }
    }
    true
}

/// coverage_of_paradox: method you already have that checks
/// each record in `enum_recs` for coverage under star logic
/// (record_matches_pop_with_stars).
pub fn coverage_of_paradox(
    sp: &Quartet<i64, i64, usize, usize>,
    enum_recs: &Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
) -> Vec<usize> {
    let mut covered = Vec::new();
    for (row_idx, rec_codes) in enum_recs.iter().enumerate() {
        let in_s1 = record_matches_pop_with_stars(rec_codes, sp.value0, size, mover);
        // let in_s1 = record_matches_pop_with_stars(rec_codes, s1, size, mover);
        let in_s2 = record_matches_pop_with_stars(rec_codes, sp.value1, size, mover);
        // let in_s2 = record_matches_pop_with_stars(rec_codes, s2, size, mover);
        if in_s1 || in_s2 {
            covered.push(row_idx);
        }
    }
    covered
}

pub fn first_equals_union_of_rest(sets: &Vec<HashSet<usize>>) -> bool {
    if sets.len() < 2 {
        // If there's only one set or none, there's nothing to compare
        return true;
    }
    
    // Create a union of all sets except the first one
    let mut union_of_rest: HashSet<usize> = HashSet::new();
    for i in 1..sets.len() {
        for &item in &sets[i] {
            union_of_rest.insert(item);
        }
    }
    
    // Check if the first set equals this union
    &sets[0] == &union_of_rest
}

pub fn coverage_of_patitioned_paradox(
    sp: &Quartet<i64, i64, usize, usize>,
    enum_recs: &Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
    lengths: &Vec<i32>,         // how many distinct codes per attribute
) -> Vec<HashSet<usize>> {
    let mut covered: Vec<HashSet<usize>> = Vec::new();

    let mut sib: HashSet<usize> = HashSet::new();
    for (row_idx, rec_codes) in enum_recs.iter().enumerate() {
        let in_s1 = record_matches_pop_with_stars(rec_codes, sp.value0, size, mover);
        let in_s2 = record_matches_pop_with_stars(rec_codes, sp.value1, size, mover);
        if in_s1 || in_s2 {
            sib.insert(row_idx);
        }
    }
    covered.push(sib);

    // For each attribute, we check if it is a star attribute
    // and if so, we partition the set of covered records
    for j in 1..lengths[sp.value2] {
        let mut sibx: HashSet<usize> = HashSet::new();
        for (row_idx, rec_codes) in enum_recs.iter().enumerate() {
            let in_s1 = record_matches_pop_with_stars(rec_codes, sp.value0 + ((j as i64) << mover[sp.value2]), size, mover);
            let in_s2 = record_matches_pop_with_stars(rec_codes, sp.value1 + ((j as i64) << mover[sp.value2]), size, mover);
            if in_s1 || in_s2 {
                sibx.insert(row_idx);
            }
        }
        covered.push(sibx);
    }

    // let result = first_equals_union_of_rest(&covered);
    // println!("Result: {}", result);

    covered
}

fn calculate_hash_for_set(set: &HashSet<usize>) -> u64 {
    let mut hasher = DefaultHasher::new();
    
    // Convert to sorted Vec to ensure consistent ordering
    let mut elements: Vec<&usize> = set.iter().collect();
    elements.sort();
    
    // Hash the sorted elements
    for element in elements {
        element.hash(&mut hasher);
    }
    
    hasher.finish()
}

pub fn coverage_sets_partitioned_equal(
    sp1: &Quartet<i64, i64, usize, usize>,
    sp2: &Quartet<i64, i64, usize, usize>,
    enum_recs: &Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
    lengths: &Vec<i32>,
) -> bool {
    let cov1 = coverage_of_patitioned_paradox(sp1, enum_recs, size, mover, lengths);
    let cov2 = coverage_of_patitioned_paradox(sp2, enum_recs, size, mover, lengths);

    // // Convert to sets and compare
    // let mut set1: HashSet<u64> = HashSet::new();
    // for hashset in &cov1 {
    //     let hash = calculate_hash_for_set(hashset);
    //     set1.insert(hash);
    // }

    // let mut set2: HashSet<u64> = HashSet::new();
    // for hashset in &cov2 {
    //     let hash = calculate_hash_for_set(hashset);
    //     set2.insert(hash);
    // }

    // Convert each HashSet to a BTreeSet, then to a HashSet of BTreeSets
    let set1: HashSet<BTreeSet<usize>> = cov1.iter()
        .map(|hash_set| hash_set.iter().cloned().collect::<BTreeSet<usize>>())
        .collect();
    
    let set2: HashSet<BTreeSet<usize>> = cov2.iter()
        .map(|hash_set| hash_set.iter().cloned().collect::<BTreeSet<usize>>())
        .collect();
    
    set1 == set2
}

/// This new method compares the coverage sets for two Simpson paradoxes
/// in the given (possibly perturbed) `enum_recs`. It returns true if
/// both paradoxes cover exactly the same records, false otherwise.
pub fn coverage_sets_equal(
    sp1: &Quartet<i64, i64, usize, usize>,
    sp2: &Quartet<i64, i64, usize, usize>,
    enum_recs: &Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
) -> bool {
    let cov1 = coverage_of_paradox(sp1, enum_recs, size, mover);
    let cov2 = coverage_of_paradox(sp2, enum_recs, size, mover);

    // Convert to sets and compare
    let set1: HashSet<usize> = cov1.into_iter().collect();
    let set2: HashSet<usize> = cov2.into_iter().collect();

    set1 == set2
}

/// coverage_of_paradox_with_perturb:
/// For each record in `enum_recs`, we check if it is covered by `sp` (star logic).
/// If covered, we add that row index to `covered`.
/// Then, with 5% probability, we modify *all* attributes in `shared_star`.
/// For each such attribute, we pick a new code in [1..lengths[attr_idx]] 
/// that is not the old value, and store it in `enum_recs[row_idx][attr_idx]`.
///
/// Returns the list of covered row indices.
pub fn coverage_of_paradox_with_perturb_sibling(
    sp: &Quartet<i64, i64, usize, usize>,
    enum_recs: &mut Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
    shared_star: &Vec<usize>,   // columns that are star attributes
    lengths: &Vec<i32>,         // how many distinct codes per attribute
){
    // let mut covered = Vec::new();
    let mut rng = rand::thread_rng();

    for (row_idx, rec) in enum_recs.iter_mut().enumerate() {
        // 1) Check coverage: record belongs to sp.value0 or sp.value1 under star logic
        let in_s1 = record_matches_pop_with_stars(rec, sp.value0, size, mover);
        let in_s2 = record_matches_pop_with_stars(rec, sp.value1, size, mover);

        if in_s1 || in_s2 {
            // The record is covered by this paradox
            // covered.push(row_idx);

            // 2) 5% probability to mutate all star attributes in this record
            if rng.gen_bool(0.05) && !shared_star.is_empty() {
                // For each star attribute in shared_star, pick a new code in [1..=max_code],
                // ensuring it is != old_val
                for &attr_idx in shared_star.iter() {
                    let max_code = lengths[attr_idx];
                    if max_code <= 1 {
                        // If max_code <= 1, we have no real alternative code
                        // or only code=1 which might be the old val
                        continue;
                    }

                    let old_val = rec[attr_idx];
                    // pick a new code != old_val
                    if old_val < 1 || old_val > max_code {
                        // old_val out of domain, pick any valid code from 1..=max_code
                        rec[attr_idx] = rng.gen_range(1..=max_code);
                    } else {
                        // old_val is in [1..max_code], pick until we differ
                        loop {
                            let candidate = rng.gen_range(1..=max_code);
                            if candidate != old_val {
                                rec[attr_idx] = candidate;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    // covered
}

pub fn coverage_of_paradox_with_perturb_division(
    sp: &Quartet<i64, i64, usize, usize>,
    enum_recs: &mut Vec<Vec<i32>>,
    size: usize,
    mover: &Vec<i32>,
    shared_star: &Vec<usize>,   // columns that are star attributes
    lengths: &Vec<i32>,         // how many distinct codes per attribute
){
    // let mut covered = Vec::new();
    let mut rng = rand::thread_rng();

    for (row_idx, rec) in enum_recs.iter_mut().enumerate() {
        // 1) Check coverage: record belongs to sp.value0 or sp.value1 under star logic
        let in_s1 = record_matches_pop_with_stars(rec, sp.value0, size, mover);
        let in_s2 = record_matches_pop_with_stars(rec, sp.value1, size, mover);

        if in_s1 || in_s2 {
            // The record is covered by this paradox
            // covered.push(row_idx);

            // 2) 5% probability to mutate all star attributes in this record
            if rng.gen_bool(0.4) && !shared_star.is_empty() {
                // For each star attribute in shared_star, pick a new code in [1..=max_code],
                // ensuring it is != old_val
                // 1) Randomly pick one attribute index from shared_star
                let attr_idx = shared_star[rng.gen_range(0..shared_star.len())];
                // for &attr_idx in shared_star.iter() {
                let max_code = lengths[attr_idx];
                if max_code <= 1 {
                    // If max_code <= 1, we have no real alternative code
                    // or only code=1 which might be the old val
                    continue;
                }

                let old_val = rec[attr_idx];
                // pick a new code != old_val
                if old_val < 1 || old_val > max_code {
                    // old_val out of domain, pick any valid code from 1..=max_code
                    rec[attr_idx] = rng.gen_range(1..=max_code);
                } else {
                    // old_val is in [1..max_code], pick until we differ
                    loop {
                        let candidate = rng.gen_range(1..=max_code);
                        if candidate != old_val {
                            rec[attr_idx] = candidate;
                            break;
                        }
                    }
                }
                // }
            }
        }
    }
    // covered
}


pub fn significance_test_driver(
    infos: &HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>>,
    enum_recs: &Vec<Vec<i32>>,  // note: we take an immutable reference now
    size: usize,
    mover: &Vec<i32>,
    lengths: &Vec<i32>,
) -> (f64, f64) {

    // We collect partial results in parallel, each partial result is (significant_count, non_significant_count).
    let results: Vec<(usize, usize)> = infos
        .par_iter()
        .map(|(group_id, group)| {
            // 1) If group <= 1 paradox, skip => return (0,0)
            if group.len() <= 1 {
                return (0usize, 0usize);
            }

            // let unique_vals: HashSet<i64> = group.iter().map(|q| q.value0).collect();
            // if unique_vals.len() <= 1 {
            //     return (0usize, 0usize);
            // }

            let unique_vals: HashSet<usize> = group.iter().map(|q| q.value2).collect();
            if unique_vals.len() <= 1 {
                return (0usize, 0usize);
            }

            // 2) Find shared_star, min_sp, max_sp
            // let (shared_star, min_sp_opt, max_sp_opt) = find_shared_star_attribute_sibling(group, size, mover);
            let (shared_star, min_sp_opt, max_sp_opt) = find_shared_star_attribute_division(group, size, mover);
            if min_sp_opt.is_none() || max_sp_opt.is_none() || group.is_empty() {
                // group invalid => (0,0)
                return (0, 0);
            }

            let sp_min = min_sp_opt.unwrap();
            let sp_max = max_sp_opt.unwrap();

            let num_runs = 50;
            let mut coverage_remain_count = 0usize;

            // Because we cannot mutate the same enum_recs across threads,
            // we clone it locally for each group.
            // We'll revert from original_data inside the same iteration, but that is local to the thread anyway.
            let original_data = enum_recs.clone();
            let mut local_data = original_data.clone();

            for _trial in 0..num_runs {
                // coverage_of_paradox_with_perturb_sibling(
                //     &sp_max,
                //     &mut local_data,
                //     size,
                //     mover,
                //     &shared_star,
                //     lengths
                // );
                coverage_of_paradox_with_perturb_division(
                    &sp_max,
                    &mut local_data,
                    size,
                    mover,
                    &shared_star,
                    lengths
                );
                // let same_coverage = coverage_sets_equal(&sp_min, &sp_max, &local_data, size, mover);
                let same_coverage = coverage_sets_partitioned_equal(&sp_min, &sp_max, &local_data, size, mover, lengths);
                if same_coverage {
                    coverage_remain_count += 1;
                }
                // revert local_data
                local_data = original_data.clone();
            }

            let fraction = coverage_remain_count as f64 / num_runs as f64;
            // println!("Group {}: {:.4}", group_id, fraction);
            if fraction <= 0.05 {
                (1, 0) // significant
            } else {
                (0, 1) // non-significant
            }
        })
        .collect();

    // Now we sum up all partial results
    let mut significant = 0usize;
    let mut non_significant = 0usize;
    for (sig, non_sig) in results {
        significant += sig;
        non_significant += non_sig;
    }

    let total = significant + non_significant;
    println!("Total: {}, Significant: {}, Non-significant: {}", total, significant, non_significant);
    if total == 0 {
        return (0.0, 0.0);
    }

    let significant_ratio = significant as f64 / total as f64;
    let non_significant_ratio = non_significant as f64 / total as f64;
    (significant_ratio, non_significant_ratio)
}


// ------------------------------------------------------------------
// Example main that uses the user arguments, calls the logic, and writes output.
// ------------------------------------------------------------------
fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: <SIZE> <NUM_LABELS> <FILENAME> [threshold]");
        return;
    }

    // Parse arguments
    let size = args[1].parse::<usize>().unwrap_or(0);
    let num_labels = args[2].parse::<usize>().unwrap_or(0);
    let filename = &args[3];
    // Optional threshold argument
    let mut threshold = if args.len() > 4 {
        args[4].parse::<f64>().unwrap_or(0.0)
    } else {
        0.0
    };

    // Start timing
    let start_time = Instant::now();

    // 1. Read CSV data
    // println!("Dataset: {}", filename);
    let lines = get_lines(filename, size);  
    // println!("Number of records: {}", lines.len());

    threshold = threshold * lines.len() as f64 / 100.0;

    // 2. Create index mappings
    let mut index_holder: HashMap<String, i32> = HashMap::new();
    let mut lengths = vec![0i32; size];             // track unique counts
    let mut enum_recs = Vec::with_capacity(lines.len());
    let mut b_labs = Vec::with_capacity(lines.len());

    create_index(
        &lines,
        size,
        num_labels,
        &mut index_holder,
        &mut lengths,
        &mut enum_recs,
        &mut b_labs,
    );

    // 3. Compute bit shifts
    let mov = get_mover(&lengths, size);

    // 4. Aggregate
    let mut aggregations: HashMap<i64, Vec<i64>> = HashMap::new();
    let mut covers: HashMap<i64, i64> = HashMap::new();

    aggregate(
        &mov,
        &lengths,
        &index_holder,
        &mut aggregations,
        &enum_recs,
        &b_labs,
        &mut covers,
        size,
        num_labels,
        threshold,
    );

    // 5. Group covers
    let grouped_covers = group_covers(&covers);
    // println!("Total number of non-empty populations: {}", covers.len());
    // println!("Total number of coverage equivalent subsets: {}", grouped_covers.len());
    // println!("Total time for materialization: {:.4} seconds", start_time.elapsed().as_secs_f64());

    // 6. Find SP
    let tup = find_sp(
        &aggregations,
        &grouped_covers,
        &lengths,
        &mov,
        &covers,
        size,
        num_labels,
    );
    let sps = tup.value0;
    let infos = tup.value1;

    // Count the number of HashSets with size >= 2
    let red = infos.iter().filter(|(_, set)| set.len() >= 2).count();
    println!("Total number of redundant relations: {}", red);
    println!("Total runtime: {:.4} seconds", start_time.elapsed().as_secs_f64());

    let count_sibling = infos.iter()
    .filter(|(_, group)| {
        // Collect .value0 from each Quartet into a set
        let unique_first_vals: std::collections::HashSet<i64> = group
            .iter()
            .map(|q| q.value0)
            .collect();
        // If that set has >=2, we say it meets our criterion
        unique_first_vals.len() >= 2
    })
    .count();

    println!("Number of sibling equivalences: {}", count_sibling);

    let count_division = infos.iter()
    .filter(|(_, group)| {
        // Collect .value0 from each Quartet into a set
        let unique_third_vals: std::collections::HashSet<usize> = group
            .iter()
            .map(|q| q.value2)
            .collect();
        // If that set has >=2, we say it meets our criterion
        unique_third_vals.len() >= 2
    })
    .count();

    println!("Number of division equivalences: {}", count_division);

    let num_runs = 100;            // how many times to repeat the random swaps
    let coverage_percentage = 0.05; // 5% of coverage for determining n_swaps

    // 7. Significance test
    let (sig_ratio, non_sig_ratio) = significance_test_driver(
        &infos,
        &mut enum_recs,
        size,
        &mov,
        &lengths,
    );

    // Convert to percentage and keep two decimal places
    println!(
        "Significant: {:.2}%, Non-significant: {:.2}%",
        sig_ratio * 100.0,
        non_sig_ratio * 100.0
    );
}