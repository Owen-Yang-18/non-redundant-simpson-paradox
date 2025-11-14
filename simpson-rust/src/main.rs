use std::io::{BufRead, BufReader};
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::Write;
// use std::path::Path;
// use std::path::PathBuf;
use std::time::Instant;
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::Entry;
use std::num::Wrapping;

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


/// Extract header from CSV file
pub fn extract_header(filename: &str) -> Vec<String> {
    if let Ok(file) = File::open(filename) {
        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        if reader.read_line(&mut first_line).is_ok() {
            return first_line
                .trim()
                .split(',')
                .map(|s| s.to_string())
                .collect();
        }
    }
    Vec::new()
}


pub fn create_index(
    records: &Vec<Vec<String>>,
    size: usize,
    num_labels: usize,
    index_holder: &mut HashMap<String, i32>,
    lengths: &mut Vec<i32>,
    enum_recs: &mut Vec<Vec<i32>>,
    binary_labels: &mut Vec<Vec<i32>>,
    reverse_index: &mut HashMap<(usize, i32), String>,
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

                    reverse_index.insert((i, lengths[i]), val.clone());
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

//     println!("Agg size: {}", agg.len());

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


/// Decode a population key to human-readable string
pub fn decode_population(
    pop_key: i64,
    mover: &Vec<i32>,
    size: usize,
    reverse_index: &HashMap<(usize, i32), String>,
    attr_names: &[String],
) -> String {
    let indices = get_indexes(pop_key, mover, size);
    
    let mut parts = Vec::new();
    for (attr_idx, &value_id) in indices.iter().enumerate() {
        let attr_name = if attr_idx < attr_names.len() {
            &attr_names[attr_idx]
        } else {
            "?"
        };
        
        if value_id == 0 {
            // Wildcard
            parts.push(format!("{}=*", attr_name));
        } else {
            // Lookup original string value
            let value_str = reverse_index
                .get(&(attr_idx, value_id as i32))
                .map(|s| s.as_str())
                .unwrap_or("?");
            parts.push(format!("{}={}", attr_name, value_str));
        }
    }
    
    format!("({})", parts.join(", "))
}


/// Write infos to human-readable file
pub fn write_infos_to_file(
    infos: &HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>>,
    filename: &str,
    mover: &Vec<i32>,
    size: usize,
    reverse_index: &HashMap<(usize, i32), String>,
    attr_names: &[String],
    label_names: &[String],
) -> std::io::Result<()> {
    use std::io::Write;
    
    let mut file = File::create(filename)?;
    
    // Sort by size (descending) and take top 100k
    let mut groups: Vec<_> = infos.iter().collect();
    groups.sort_by(|a, b| {
        // Sort by set size (descending), then by signature for consistency
        b.1.len().cmp(&a.1.len())
            .then_with(|| a.0.cmp(b.0))
    });
    
    let max_groups = 100_000;
    let groups_to_write = groups.len().min(max_groups);
    let omitted = groups.len().saturating_sub(max_groups);
    
    writeln!(file, "Total redundant paradox groups: {}", infos.len())?;
    writeln!(file, "Showing top {} groups (by size)", groups_to_write)?;
    if omitted > 0 {
        writeln!(file, "Omitted {} smaller groups", omitted)?;
    }
    writeln!(file, "{}", "=".repeat(80))?;
    writeln!(file)?;
    
    for (group_idx, (signature, paradox_set)) in groups.iter().take(max_groups).enumerate() {
        writeln!(file, "Group {} (Signature: {}, Count: {})", 
                 group_idx + 1, signature, paradox_set.len())?;
        writeln!(file, "{}", "-".repeat(80))?;
        
        let mut paradoxes: Vec<_> = paradox_set.iter().collect();
        paradoxes.sort_by_key(|sp| (sp.value0, sp.value1, sp.value2, sp.value3));
        
        for (sp_idx, sp) in paradoxes.iter().enumerate() {
            // Decode populations
            let s1_str = decode_population(sp.value0, mover, size, reverse_index, attr_names);
            let s2_str = decode_population(sp.value1, mover, size, reverse_index, attr_names);
            
            // Get separator name
            let sep_name = if sp.value2 < attr_names.len() {
                &attr_names[sp.value2]
            } else {
                "?"
            };
            
            // Get label name
            let label_name = if sp.value3 < label_names.len() {
                &label_names[sp.value3]
            } else {
                "?"
            };
            
            writeln!(file, "  Paradox {}: s1={}, s2={}, separator={}, label={}",
                     sp_idx + 1, s1_str, s2_str, sep_name, label_name)?;
        }
        writeln!(file)?;
    }
    
    Ok(())
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
    let threshold = if args.len() > 4 {
        args[4].parse::<f64>().unwrap_or(0.0)
    } else {
        0.0
    };

    // Start timing
    let start_time = Instant::now();

    // 1. Read CSV data
    println!("Dataset: {}", filename);
    let lines = get_lines(filename, size);  
    println!("Number of records: {}", lines.len());

    let header = extract_header(filename);
    let attr_names: Vec<String> = if header.len() >= size {
        header[0..size].to_vec()
    } else {
        (0..size).map(|i| format!("attr_{}", i)).collect()
    };
    let label_names: Vec<String> = if header.len() > size {
        header[size..].to_vec()
    } else {
        (0..num_labels).map(|i| format!("label_{}", i)).collect()
    };

    // 2. Create index mappings
    let mut index_holder: HashMap<String, i32> = HashMap::new();
    let mut lengths = vec![0i32; size];             // track unique counts
    let mut enum_recs = Vec::with_capacity(lines.len());
    let mut b_labs = Vec::with_capacity(lines.len());

    let mut reverse_index: HashMap<(usize, i32), String> = HashMap::new();

    create_index(
        &lines,
        size,
        num_labels,
        &mut index_holder,
        &mut lengths,
        &mut enum_recs,
        &mut b_labs,
        &mut reverse_index,
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
    println!("Total time for materialization: {:.4} seconds", start_time.elapsed().as_secs_f64());
    println!("Total number of non-empty populations: {}", covers.len());
    println!("Total number of coverage equivalent subsets: {}", grouped_covers.len());

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

    println!("Total runtime: {:.4} seconds", start_time.elapsed().as_secs_f64());
    println!("Total number of SP: {}", sps.len());
    println!("Total number of redundant relations: {}", infos.len());

    use std::path::Path;
    
    let output_filename = {
        let input_path = Path::new(filename);
        let stem = input_path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        format!("{}_infos.txt", stem)
    };
    
    match write_infos_to_file(
        &infos,
        &output_filename,
        &mov,
        size,
        &reverse_index,
        &attr_names,
        &label_names,
    ) {
        Ok(_) => println!("Results written to {}", output_filename),
        Err(e) => eprintln!("Failed to write output file: {}", e),
    }
}
