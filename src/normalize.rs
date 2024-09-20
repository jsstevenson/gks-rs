enum Mode {
    Expand,
    LeftShuffle,
    RightShuffle,
    TrimOnly,
    Vcf,
}

fn trim_left<'a>(alleles: Vec<&'a str>) -> (usize, Vec<&'a str>) {
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

    let new_alleles: Vec<&'a str> = alleles.iter().map(|x| &x[trimmed..]).collect();
    (trimmed, new_alleles)
}

fn trim_right<'a>(alleles: Vec<&'a str>) -> (usize, Vec<&'a str>) {
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

    let new_alleles: Vec<&'a str> = alleles
        .iter()
        .map(|x| &x[..(x.len() - trimmed)])
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

fn normalize<'a>(
    sequence: &str,
    mut interval: (usize, usize),
    alleles: Vec<&'a str>,
    mut bounds: Option<(usize, usize)>,
    mut mode: Option<Mode>,
    mut anchor_length: Option<usize>,
    mut trim: Option<bool>,
) -> Result<((usize, usize), Vec<&'a str>), Box<dyn std::error::Error>> {
    if let Some(_) = bounds {
    } else {
        bounds = Some((0, sequence.len()))
    };
    if let Some(_) = mode {
    } else {
        mode = Some(Mode::Expand)
    };
    if let Some(_) = anchor_length {
    } else {
        anchor_length = Some(0)
    };
    if let Some(_) = trim {
    } else {
        trim = Some(true)
    };
    let (start, end) = interval;
    if start > end {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "Invalid interval: end ({}) is less than start ({})",
                end, start
            ),
        )));
    }

    let left_anchor;
    let right_anchor;
    if let Some(Mode::Vcf) = mode {
        if anchor_length.unwrap() == 0 {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("May not provide non-zero anchor size with VCF normalization mode"),
            )));
        }
        if !trim.unwrap()  {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("May not disable trimming with VCF normalization mode"),
            )));
        }
        left_anchor = 1;
        right_anchor = 0;
        mode = Some(Mode::LeftShuffle);
    } else {
        left_anchor = anchor_length.unwrap();
        right_anchor = anchor_length.unwrap();
    }

    // TODO raises value error here bc of requirement that first allele (REF) is None
    // seems unnecessary...

    if trim.unwrap() {
        let (left_trimmed, alleles) = trim_left(alleles);
        interval.0 += left_trimmed;
        let (right_trimmed, alleles) = trim_right(alleles);
        interval.1 -= right_trimmed;
    }

    let lens: Vec<usize> = alleles.iter().map(|a| a.len()).collect();


    Ok(((0, 1), vec![""]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize() {
        let sequence = "CCCCCCCCACACACACACTAGCAGCAGCA";
        let (bounds, normalized_alleles) = normalize(
            sequence,
            (22, 25),
            vec!["", "GC", "AGCAC"],
            None,
            Some(Mode::TrimOnly),
            None,
            None,
        )
        .unwrap();
        assert_eq!(bounds, (22, 24));
        assert_eq!(normalized_alleles, vec!["AG", "G", "AGCA"]);

        let (bounds, normalized_alleles) = normalize(
            sequence,
            (22, 22),
            vec!["", "AGC"],
            None,
            Some(Mode::RightShuffle),
            None,
            None,
        )
        .unwrap();
        assert_eq!(bounds, (29, 29));
        assert_eq!(normalized_alleles, vec!["", "GCA"]);

        let (bounds, normalized_alleles) = normalize(
            sequence,
            (22, 22),
            vec!["", "AGC"],
            None,
            Some(Mode::Expand),
            None,
            None,
        )
        .unwrap();
        assert_eq!(bounds, (19, 29));
        assert_eq!(normalized_alleles, vec!["AGCAGCAGCA", "AGCAGCAGCAGCA"]);
    }

    #[test]
    fn test_trim_left() {
        let mut alleles: Vec<&str> = vec!["", "AA"];
        let mut num_trimmed: usize;
        let mut new_alleles: Vec<&str>;
        (num_trimmed, new_alleles) = trim_left(alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["".to_string(), "AA".to_string()]);

        alleles = vec!["A", "AA"];
        (num_trimmed, new_alleles) = trim_left(alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["".to_string(), "A".to_string()]);

        alleles = vec!["AT", "AA"];
        (num_trimmed, new_alleles) = trim_left(alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["T".to_string(), "A".to_string()]);

        alleles = vec!["AA", "AA"];
        (num_trimmed, new_alleles) = trim_left(alleles);
        assert_eq!(num_trimmed, 2);
        assert_eq!(new_alleles, vec!["".to_string(), "".to_string()]);

        alleles = vec!["CAG", "CG"];
        (num_trimmed, new_alleles) = trim_left(alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["AG".to_string(), "G".to_string()]);
    }

    #[test]
    fn test_trim_right() {
        let mut alleles: Vec<&str> = vec!["", "AA"];
        let mut num_trimmed: usize;
        let mut new_alleles: Vec<&str>;
        (num_trimmed, new_alleles) = trim_right(alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["".to_string(), "AA".to_string()]);

        alleles = vec!["A", "AA"];
        (num_trimmed, new_alleles) = trim_right(alleles);
        assert_eq!(num_trimmed, 1);
        assert_eq!(new_alleles, vec!["".to_string(), "A".to_string()]);

        alleles = vec!["AT", "AA"];
        (num_trimmed, new_alleles) = trim_right(alleles);
        assert_eq!(num_trimmed, 0);
        assert_eq!(new_alleles, vec!["AT".to_string(), "AA".to_string()]);

        alleles = vec!["AA", "AA"];
        (num_trimmed, new_alleles) = trim_right(alleles);
        assert_eq!(num_trimmed, 2);
        assert_eq!(new_alleles, vec!["".to_string(), "".to_string()]);

        alleles = vec!["CAG", "CG"];
        (num_trimmed, new_alleles) = trim_right(alleles);
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
