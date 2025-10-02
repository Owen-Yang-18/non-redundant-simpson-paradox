use std::io::{BufRead, BufReader};
// use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::time::Instant;
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::Entry;
use std::num::Wrapping;
use std::cmp::{max};
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


pub fn aggregate(
    mov: &Vec<i32>,
    _lengths: &Vec<i32>,
    _index_holder: &HashMap<String, i32>, // not deeply used in your snippet
    map: &mut HashMap<i64, Vec<i64>>,
    enum_recs: &Vec<Vec<i32>>,
    binary_labels: &Vec<Vec<i32>>,
    covers: &mut HashMap<i64, i64>,
    size: usize,
    num_labels: usize,
    threshold: f64,
) {
    let mut agg: HashMap<Vec<i32>, Vec<i64>> = HashMap::new();
    let hash: i64 = 1;
    let mut count = 0i64;
    let large_prime: i64 = 2305843009213693951; // 2^61 - 1

    let mut bases: HashMap<Vec<i32>, i64> = HashMap::new();

    for list_key in enum_recs {
        if !agg.contains_key(list_key) {
            // create new label array
            let mut labels = vec![0i64; num_labels];
            for y in 0..num_labels {
                // in Java: (1L << 32) + binaryLabels[count][y]
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
            // update existing label array
            let existing = agg.get_mut(list_key).unwrap();
            for y in 0..num_labels {
                existing[y] += (1i64 << 32) + (binary_labels[count as usize][y] as i64);
            }
        }
        count += 1;
    }

    // println!("Agg size: {}", agg.len());

    for i in 0..=size {
        let mut cov: HashMap<i64, i64> = HashMap::new();
        let mut stat: HashMap<i64, Vec<i64>> = HashMap::new();

        for list_key in agg.keys() {
            let base = bases[list_key];
            let labels = &agg[list_key];

            if i == 0 {
                // everything merges into key 0
                let l = 0i64;
                map.entry(l).or_insert(vec![0; num_labels]);
                covers.entry(l).or_insert(hash);
                
                let entry_map = map.get_mut(&l).unwrap();
                for y in 0..num_labels {
                    entry_map[y] += labels[y];
                }

                let cv = covers.get_mut(&l).unwrap();
                *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64; 

            } else if i == size {
                // only add if #inc >= threshold
                // from Java: (stat.get(l)[0] >> 32) is the total
                // (stat.get(l)[0] % (1 << 32)) is inc but the code is a bit different
                // Actually: getRate(long) => inc / tot
                // We check if tot >= threshold or inc >= threshold?
                // The Java snippet does: if (double)(agg.get(listKey)[0] >> 32) >= threshold
                // That means if total >= threshold. 
                let tot = (labels[0] >> 32) as f64;
                if tot >= threshold {
                    map.insert(base, labels.clone());
                    let cval = (hash * (base + 31)) % large_prime;
                    covers.insert(base, cval);
                }

            } else {
                // loop over subsets of size i
                // mask from (1 << i) - 1 up to ...
                let start_mask: i64 = ((1i64 << i) - 1) << (size - i);
                // The Java code: for (long mask = (1L << i) - 1; mask <= ((1L << i) - 1) << (SIZE - i); ...)
                // We must be careful about how we iterate in Rust because we do next_mask.
                // So we do something like:

                let mut mask = (1i64 << i) - 1;
                while mask <= start_mask {
                    let mut l = 0i64;
                    for j in 0..size {
                        let ind = 1i64 << j;
                        if (ind & mask) == ind {
                            l += (list_key[j] as i64) << mov[j];
                        }
                    }

                    // now update stat and cov
                    match stat.entry(l) {
                        Entry::Vacant(e) => {
                            e.insert(labels.clone());
                            cov.insert(l, hash);
                            // update covers
                            let cv = cov.get_mut(&l).unwrap();
                            // *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64;
                            *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

                            // check threshold
                            let tot = (labels[0] >> 32) as f64;
                            if tot >= threshold {
                                map.insert(l, labels.clone());
                                covers.insert(l, *cv);
                            }
                        }
                        Entry::Occupied(mut e) => {
                            // If already in map but not in the global map
                            // or partially in the global map
                            let e_val = e.get_mut();
                            for y in 0..num_labels {
                                e_val[y] += labels[y];
                            }

                            let cv = cov.get_mut(&l).unwrap();
                            // *cv = ((*cv as i128) * ((base + 31) as i128) % (large_prime as i128)) as i64;
                            *cv = (Wrapping(*cv) * Wrapping(base + 31)).0 % large_prime;

                            let tot = (e_val[0] >> 32) as f64;
                            if tot >= threshold {
                                // update global
                                map.insert(l, e_val.clone());
                                covers.insert(l, *cv);
                            }
                        }
                    }

                    // proceed to next mask
                    mask = next_mask(mask);
                }
            }
        }
    }
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


pub fn find_sp(
    aggregations: &HashMap<i64, Vec<i64>>,
    grouped_covers: &HashMap<i64, Vec<i64>>,   // changed here
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
    let mut all_sps: HashSet<Quartet<i64, i64, usize, usize>> = HashSet::new();
    let mut infos: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();
    let mut infos_aux: HashMap<i64, HashSet<Quartet<i64, i64, usize, usize>>> = HashMap::new();

    // groupedCovers: key = cover value, value = Vec<i64> of populations
    for (_cover_val, pop_list) in grouped_covers {
        // getLeastStars now uses &Vec<i64>
        let stars = get_least_stars(pop_list, mover, size);
        let l1 = stars.value1;
        let zeros = stars.value0;

        if zeros.len() < 2 {
            continue;
        }

        // for each attribute i in zeros
        for &i in zeros.iter() {
            // i1 in [1..lengths[i]]
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
                            if all_sps.contains(&p) {
                                continue;
                            }
                            let res = check_sp(c1, c2, j, y, aggregations, mover, lengths, covers);
                            let is_sp = res.value0;
                            let info = res.value1;

                            if is_sp && !infos_aux.contains_key(&info) {
                                // create a new "info bucket"
                                all_sps.insert(p.clone());
                                let mut new_set = HashSet::new();
                                new_set.insert(p.clone());
                                infos_aux.insert(info.clone(), new_set.clone());
                                infos.insert(info.clone(), new_set);

                                // find siblings
                                let c1cov = covers.get(&c1).unwrap();
                                let c2cov = covers.get(&c2).unwrap();

                                let c1cov_list = grouped_covers.get(c1cov).unwrap();
                                let c2cov_list = grouped_covers.get(c2cov).unwrap();

                                // examine all pairs in c1cov_list x c2cov_list
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
                                            let big = if t1 < t2 { t2 } else { t1 };

                                            let p_prime = Quartet::new(small, big, j, y);
                                            all_sps.insert(p_prime.clone());

                                            let info_set = infos.get_mut(&info).unwrap();
                                            info_set.insert(p_prime.clone());
                                            let info_aux_set = infos_aux.get_mut(&info).unwrap();
                                            info_aux_set.insert(p_prime);
                                        }
                                    }
                                }
                            } else if is_sp {
                                // If is_sp and info already present, copy in the new p or merge
                                if let Some(old_set) = infos_aux.get(&info) {
                                    for p_prime in old_set.iter() {
                                        let p_hat = Quartet::new(p_prime.value0, p_prime.value1, j, y);
                                        all_sps.insert(p_hat.clone());
                                        if let Some(infos_set) = infos.get_mut(&info) {
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
    }

    Pair::new(all_sps, infos)
}


fn decode_counts(val: i64) -> (i64, i64) {
    let tot = val >> 32;
    let inc = val & 0xFFFF_FFFF;
    (inc, tot)
}

fn encode_counts(inc: i64, tot: i64) -> i64 {
    (tot << 32) | (inc & 0xFFFF_FFFF)
}

fn sum_coverage_for_paradox(
    s1: i64,
    s2: i64,
    j: usize,
    y: usize,
    aggregations: &HashMap<i64, Vec<i64>>,
    lengths: &Vec<i32>,
    mov: &Vec<i32>,
) -> i64 {
    let max_k = lengths[j] as i64;
    let mut total = 0i64;
    for k in 1..=max_k {
        let s1_j = s1 + (k << mov[j]);
        let s2_j = s2 + (k << mov[j]);

        if let Some(vec_s1) = aggregations.get(&s1_j) {
            let (_, tot1) = decode_counts(vec_s1[y]);
            total += tot1;
        }
        if let Some(vec_s2) = aggregations.get(&s2_j) {
            let (_, tot2) = decode_counts(vec_s2[y]);
            total += tot2;
        }
    }
    total
}


/// Randomly transfers one label=1 from s1_j -> s2_j or from s2_j -> s1_j,
/// if constraints allow (inc1>0, inc2<tot2, etc.). Uses separate scopes
/// to avoid multiple mutable borrows on aggregator.
// fn transfer_one_label(
//     aggregator: &mut HashMap<i64, Vec<i64>>,
//     s1_j: i64,
//     s2_j: i64,
//     y: usize,
// ) {
//     // Step 1: read inc1, tot1 from s1_j in its own scope
//     let mut inc1 = 0;
//     let mut tot1 = 0;
//     {
//         if let Some(s1_vec) = aggregator.get_mut(&s1_j) {
//             let (i, t) = decode_counts(s1_vec[y]);
//             inc1 = i;
//             tot1 = t;
//         }
//     }

//     // Step 2: read inc2, tot2 from s2_j in its own scope
//     let mut inc2 = 0;
//     let mut tot2 = 0;
//     {
//         if let Some(s2_vec) = aggregator.get_mut(&s2_j) {
//             let (i, t) = decode_counts(s2_vec[y]);
//             inc2 = i;
//             tot2 = t;
//         }
//     }

//     // Step 3: decide transfer direction via coin flip
//     // If direction=0 => s1_j -> s2_j, else s2_j -> s1_j
//     let mut rng = rand::thread_rng();
//     let direction = rng.gen_bool(0.5);

//     if direction {
//         // Attempt to move one label from s1_j to s2_j
//         if inc1 > 0 && inc2 < tot2 {
//             inc1 -= 1;
//             inc2 += 1;
//         }
//     } else {
//         // Attempt to move one label from s2_j to s1_j
//         if inc2 > 0 && inc1 < tot1 {
//             inc2 -= 1;
//             inc1 += 1;
//         }
//     }

//     // Step 4: write back inc1, tot1
//     {
//         if let Some(s1_vec) = aggregator.get_mut(&s1_j) {
//             s1_vec[y] = encode_counts(inc1, tot1);
//         }
//     }

//     // Step 5: write back inc2, tot2
//     {
//         if let Some(s2_vec) = aggregator.get_mut(&s2_j) {
//             s2_vec[y] = encode_counts(inc2, tot2);
//         }
//     }
// }

fn transfer_one_label(
    aggregator: &mut HashMap<i64, Vec<i64>>,
    s1_j: i64,
    s2_j: i64,
    y: usize,
) {
    // We will read the inc/tot from s1_j, store them in local variables,
    // then in a separate scope read from s2_j, store them, and do the logic.

    // 1) Pull out s1_j (in its own scope), read inc1, tot1
    let mut inc1 = 0;
    let mut tot1 = 0;
    {
        if let Some(s1_vec) = aggregator.get_mut(&s1_j) {
            let (i, t) = decode_counts(s1_vec[y]);
            inc1 = i;
            tot1 = t;
            // If you also needed to *modify* s1_vec[y] right now, do it here.
            // e.g. s1_vec[y] = encode_counts(i+10, t+10);
        }
    } // s1_vec goes out of scope, so the mutable borrow on aggregator for s1_j is dropped here.

    // 2) Pull out s2_j (in its own scope), read inc2, tot2
    let mut inc2 = 0;
    let mut tot2 = 0;
    {
        if let Some(s2_vec) = aggregator.get_mut(&s2_j) {
            let (i, t) = decode_counts(s2_vec[y]);
            inc2 = i;
            tot2 = t;
            // Similarly, if you needed to modify s2_vec[y] immediately, do it here.
        }
    } // s2_vec goes out of scope, dropping the borrow for s2_j.

    // 3) We now have inc1, tot1, inc2, tot2 in local variables.
    //    Decide if we can "transfer one label=1" from s1_j -> s2_j.
    if inc1 > 0 && inc2 < tot2 && s1_j != s2_j {
        inc1 -= 1;  // "move" one record from s1_j
        inc2 += 1;  // "add" that record to s2_j
    }

    // 4) Write back the new inc1, tot1 to s1_j (again in a new scope).
    {
        if let Some(s1_vec) = aggregator.get_mut(&s1_j) {
            s1_vec[y] = encode_counts(inc1, tot1);
        }
    }

    // 5) Write back the new inc2, tot2 to s2_j (final scope).
    {
        if let Some(s2_vec) = aggregator.get_mut(&s2_j) {
            s2_vec[y] = encode_counts(inc2, tot2);
        }
    }
}


fn perturb_paradox_n_swaps(
    s1: i64,
    s2: i64,
    j: usize,
    y: usize,
    n: usize,
    aggregations: &HashMap<i64, Vec<i64>>,
    // check_sp: CheckSpFn,
    // covers: &HashMap<i64, i64>,
    mov: &Vec<i32>,
    lengths: &Vec<i32>,
) -> bool {
    let mut aggregator_copy = aggregations.clone();

    // All sub-pop pairs
    let max_k = lengths[j] as i64;
    let mut pairs = Vec::with_capacity(max_k as usize);
    for idx in 1..=max_k {
        let s1_j = s1 + (idx << mov[j]);
        let s2_j = s2 + (idx << mov[j]);
        pairs.push((s1_j, s2_j));
    }

    // Randomly do n single transfers
    let mut rng = rand::thread_rng();
    for _ in 0..n {
        let pick = rng.gen_range(0..pairs.len());
        let (sub1, sub2) = pairs[pick];
        transfer_one_label(&mut aggregator_copy, sub1, sub2, y);
    }

    // Now check if it remains a paradox
    let result = check_sp_sign(s1, s2, j, y, &aggregator_copy, mov, lengths);
    result
}


fn sp_significance_single(
    s1: i64,
    s2: i64,
    j: usize,
    y: usize,
    aggregations: &HashMap<i64, Vec<i64>>,
    // covers: &HashMap<i64, i64>,
    mov: &Vec<i32>,
    lengths: &Vec<i32>,
    // check_sp: CheckSpFn,
    num_runs: usize,
    coverage_percentage: f64, // e.g. 0.05 for 5%
) -> f64 {
    // 1) Calculate coverage
    let total_coverage = sum_coverage_for_paradox(s1, s2, j, y, aggregations, lengths, mov);
    // multiply by coverage_percentage, divide by 2, round up
    let n_swaps_f = (total_coverage as f64) * coverage_percentage / 2.0;
    let n_swaps = n_swaps_f.ceil() as usize;
    // let n_swaps = 2;

    // 2) Repeat the random-swaps test
    let mut remain_count = 0;
    for _ in 0..num_runs {
        let still_paradox = perturb_paradox_n_swaps(
            s1, s2, j, y, 
            n_swaps,
            aggregations,
            // check_sp,
            // covers,
            mov,
            lengths,
        );
        if still_paradox {
            remain_count += 1;
        }
    }

    // fraction of times we remain a paradox
    (remain_count as f64) / (num_runs as f64)
}


pub fn sp_significance_for_all(
    sps: &HashSet<Quartet<i64, i64, usize, usize>>,
    aggregations: &HashMap<i64, Vec<i64>>,
    // covers: &HashMap<i64, i64>,
    mov: &Vec<i32>,
    lengths: &Vec<i32>,
    // check_sp: CheckSpFn,
    num_runs: usize,
    coverage_percentage: f64,
) -> SignificanceSummary {
    // Two vectors to store coverage of significant and not-significant SPs
    let mut coverage_significant = Vec::new();
    let mut coverage_not_significant = Vec::new();

    for sp in sps.iter() {
        let s1 = sp.value0;
        let s2 = sp.value1;
        let j  = sp.value2;
        let y  = sp.value3;
        let max_k = lengths[j] as i64;

        // 1) Compute coverage
        let total_coverage = sum_coverage_for_paradox(s1, s2, j, y, aggregations, lengths, mov);

        let mut local_map = HashMap::new();

        // Insert the aggregator data for s1, s2 themselves
        if let Some(v) = aggregations.get(&s1) {
            local_map.insert(s1, v.clone());
        }
        if let Some(v) = aggregations.get(&s2) {
            local_map.insert(s2, v.clone());
        }

        for idx in 1..=max_k {
            let s1_j = s1 + (idx << mov[j]);
            let s2_j = s2 + (idx << mov[j]);

            // if aggregator has a vector for s1_j, clone it and store
            if let Some(v) = aggregations.get(&s1_j) {
                local_map.insert(s1_j, v.clone());
            }
            // likewise for s2_j
            if let Some(v) = aggregations.get(&s2_j) {
                local_map.insert(s2_j, v.clone());
            }
        }
        
        // 2) Compute significance fraction
        let fraction = sp_significance_single(
            s1, s2, j, y,
            &local_map,
            // covers,
            mov,
            lengths,
            // check_sp,
            num_runs,
            coverage_percentage
        );

        // 3) Classify
        if fraction < 0.05 {
            coverage_significant.push(total_coverage as f64);
        } else {
            coverage_not_significant.push(total_coverage as f64);
        }
    }

    // Summaries
    let num_significant = coverage_significant.len();
    let num_not_significant = coverage_not_significant.len();

    let (mean_sig, std_sig) = mean_and_std(&coverage_significant);
    let (mean_not, std_not) = mean_and_std(&coverage_not_significant);

    SignificanceSummary {
        num_significant,
        num_not_significant,

        mean_cov_significant: mean_sig,
        std_cov_significant: std_sig,

        mean_cov_not_significant: mean_not,
        std_cov_not_significant: std_not,
    }
}

/// Simple utility to compute mean & std of a slice of f64
fn mean_and_std(values: &Vec<f64>) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let n = values.len() as f64;
    let sum: f64 = values.iter().sum();
    let mean = sum / n;
    let var_sum: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum();
    let var = var_sum / n; // population variance
    let std = var.sqrt();

    (mean, std)
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

    println!("Total number of SP: {}", sps.len());
    // println!("Total number of redundant relations: {}", infos.len());
    println!("Total runtime: {:.4} seconds", start_time.elapsed().as_secs_f64());

    let num_runs = 100;            // how many times to repeat the random swaps
    let coverage_percentage = 0.1; // 10% of coverage for determining n_swaps

    let summary = sp_significance_for_all(
        &sps,
        &aggregations,
        // &covers,
        &mov,
        &lengths,
        // check_sp,   // your function pointer for check_sp
        num_runs,
        coverage_percentage
    );

    // println!("Number of significant SPs: {}", summary.num_significant);
    // println!("Mean coverage of significant: {:.2}", summary.mean_cov_significant);
    // println!("Std coverage of significant:  {:.2}", summary.std_cov_significant);

    // println!("Number of NOT significant SPs: {}", summary.num_not_significant);
    // println!("Mean coverage of not-significant: {:.2}", summary.mean_cov_not_significant);
    // println!("Std coverage of not-significant:  {:.2}", summary.std_cov_not_significant);

    println!("Percentage of signficant SPs: {:.2}%", (summary.num_significant as f64) / (sps.len() as f64) * 100.0);
    // println!("Percentage of not-significant SPs: {:.2}%", (summary.num_not_significant as f64) / (sps.len() as f64) * 100.0);
}