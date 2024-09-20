enum Mode {
    Expand,
    LeftShuffle,
    RightShuffle,
    TrimOnly,
    Vcf,
}

fn trim_left(alleles: &[String]) -> (usize, Vec<String>) {
    if alleles.is_empty() {
        return (0, vec![]);
    }

    let mut trimmed = 0;
    let min_len = alleles.iter().map(|x| x.len()).min().unwrap();

    while trimmed < min_len {
        let nexts: Vec<char> = alleles
            .iter()
            .map(|x| x.chars().nth(trimmed).unwrap())
            .collect();
        if nexts.iter().all(|&x| x == nexts[0]) {
            trimmed += 1;
        } else {
            break;
        }
    }

    let new_alleles: Vec<String> = alleles.iter().map(|x| x[trimmed..].to_string()).collect();
    (trimmed, new_alleles)
}

fn trim_right(alleles: &[String]) -> (usize, Vec<String>) {
    if alleles.is_empty() {
        return (0, vec![]);
    }

    let mut trimmed = 0;
    let min_len = alleles.iter().map(|x| x.len()).min().unwrap(); // TODO is len() inappropriate
                                                                  // here?

    while trimmed < min_len {
        let nexts: Vec<char> = alleles
            .iter()
            .map(|x| x.chars().nth(x.chars().count() - trimmed - 1).unwrap())
            .collect();
        if nexts.iter().all(|&x| x == nexts[0]) {
            trimmed += 1;
        } else {
            break;
        }
    }

    let new_alleles: Vec<String> = alleles
        .iter()
        .map(|x| x[..(x.len() - trimmed)].to_string())
        .collect();
    (trimmed, new_alleles)
}

// TODO abstract sequence interface
// probably needs more tests... feels like some consideration of bounds-based errors?
fn roll_left(sequence: &String, alleles: &[String], ref_pos: usize, bound: usize) -> usize {
    let lens: Vec<usize> = alleles.iter().map(|a| a.len()).collect();
    let mut d = 0;
    let max_d = ref_pos - bound;

    while d <= max_d
        && !alleles.iter().enumerate().any(|(i, a)| {
            let len = lens[i];
            !a.is_empty()
                && a.chars().nth((len + (len - (d + 1) % len)) % len)
                    != sequence.chars().nth(ref_pos - d)
        })
    {
        d += 1;
    }

    d
}

// TODO tests??
fn roll_right(sequence: &String, alleles: &[String], ref_pos: usize, bound: usize) -> usize {
    let max_d = bound - ref_pos;
    let lens: Vec<usize> = alleles.iter().map(|a| a.len()).collect();
    let mut d = 0;

    while d <= max_d
        && !alleles.iter().enumerate().any(|(i, a)| {
            !a.is_empty() && a.chars().nth(d % lens[i]) != sequence.chars().nth(ref_pos + d)
        })
    {
        d += 1;
    }

    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_left() {
        let mut alleles: Vec<String> = vec![String::from(""), String::from("AA")];
        let mut num_trimmed: usize;
        let mut new_alleles: Vec<String>;
        (num_trimmed, new_alleles) = trim_left(&alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["".to_string(), "AA".to_string()]);

        alleles = vec![String::from("A"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_left(&alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["".to_string(), "A".to_string()]);

        alleles = vec![String::from("AT"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_left(&alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["T".to_string(), "A".to_string()]);

        alleles = vec![String::from("AA"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_left(&alleles);
        assert_eq!(num_trimmed, 2);
        assert_eq!(new_alleles, vec!["".to_string(), "".to_string()]);

        alleles = vec![String::from("CAG"), String::from("CG")];
        (num_trimmed, new_alleles) = trim_left(&alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["AG".to_string(), "G".to_string()]);
    }

    #[test]
    fn test_trim_right() {
        let mut alleles: Vec<String> = vec![String::from(""), String::from("AA")];
        let mut num_trimmed: usize;
        let mut new_alleles: Vec<String>;
        (num_trimmed, new_alleles) = trim_right(&alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["".to_string(), "AA".to_string()]);

        alleles = vec![String::from("A"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_right(&alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["".to_string(), "A".to_string()]);

        alleles = vec![String::from("AT"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_right(&alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["AT".to_string(), "AA".to_string()]);

        alleles = vec![String::from("AA"), String::from("AA")];
        (num_trimmed, new_alleles) = trim_right(&alleles);
        assert_eq!(num_trimmed, 2);
        assert_eq!(new_alleles, vec!["".to_string(), "".to_string()]);

        alleles = vec![String::from("CAG"), String::from("CG")];
        (num_trimmed, new_alleles) = trim_right(&alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["CA".to_string(), "C".to_string()]);
    }

    #[test]
    fn test_roll_left() {
        let sequence = "ACGTACGT".to_string();
        let mut alleles: Vec<String> = vec![String::from("ACG"), String::from("ACG")];
        let mut roll_index: usize;
        roll_index = roll_left(&sequence, &alleles, 6, 0);
        assert_eq!(roll_index, 3);

        alleles = vec![String::from("ACGT")];
        roll_index = roll_left(&sequence, &alleles, 7, 0);
        assert_eq!(roll_index, 8);

        roll_index = roll_left(&sequence, &alleles, 7, 3);
        assert_eq!(roll_index, 5);
    }

    #[test]
    fn test_roll_right() {
        let sequence = "ACGTACGT".to_string();
        let mut alleles: Vec<String> = vec![String::from("ACGT"), String::from("ACGT")];
        let mut roll_index: usize;
        roll_index = roll_right(&sequence, &alleles, 0, 3);
        assert_eq!(roll_index, 4);

        roll_index = roll_right(&sequence, &alleles, 0, 7);
        assert_eq!(roll_index, 8);

        alleles = vec![String::from("ACG")];
        roll_index = roll_right(&sequence, &alleles, 0, 7);
        assert_eq!(roll_index, 3);

        alleles = vec![String::from(""), String::from("ACGT")];
        roll_index = roll_right(&sequence, &alleles, 0, 7);
        assert_eq!(roll_index, 8);

    }
}
