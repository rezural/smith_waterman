#![feature(test)]
extern crate test;
extern crate smith_waterman;
mod tests{
    use smith_waterman::*;
    use test::Bencher;
    fn genome_sequence_many(count: isize) -> String{
        let mut mill = "AACACACT".to_string();
        for _ in 1..count{
            mill.push_str("AACACACT");
        }
        mill
    }

    fn read_sequence_many(count: isize) -> String{
        let mut mill = "AAGCACAC".to_string();
        for _ in 1..count{
            mill.push_str("AAGCACAC");
        }
        mill
    }
    fn genome_sequence() -> String{
        "AACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTAACACACTACACACTA".to_string()
    }

    fn read_sequence() -> String{
        "AAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAAGCACACAGCACACA".to_string()
    }
    #[bench]
    fn bench_maxtrix_loops(b: &mut Bencher) {
    return;
        let g = genome_sequence().chars().collect();
        let r = genome_sequence().chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_loops());
    }

    #[bench]
    fn bench_maxtrix_loops_many(b: &mut Bencher) {
    return;
        let g = genome_sequence_many(1_000).chars().collect();
        let r = read_sequence_many(1_000).chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_loops());
    }
    #[bench]
    fn bench_maxtrix_loops_many_many(b: &mut Bencher) {
    return;
        let g = genome_sequence_many(10_000).chars().collect();
        let r = read_sequence_many(1_000).chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_loops());
    }
    #[bench]
    fn bench_maxtrix_fn(b: &mut Bencher) {
        let g = genome_sequence().chars().collect();
        let r = genome_sequence().chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_fn());
    }
    #[bench]
    fn bench_maxtrix_fn_many(b: &mut Bencher) {
    return;
        let g = genome_sequence_many(1_000).chars().collect();
        let r = read_sequence_many(1_000).chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_fn());
    }
    #[bench]
    fn bench_maxtrix_fn_many_many(b: &mut Bencher) {
    return;

        let g = genome_sequence_many(10_000).chars().collect();
        let r = read_sequence_many(1_000).chars().collect();
        let mut smitty = SmithWaterman::new(&g, &r);
        b.iter(|| smitty.set_matrix_fn());
    }
    #[bench]
    fn bench_maxtrix_thread(b: &mut Bencher) {
        let g: Vec<char>= genome_sequence().chars().collect();
        let r: Vec<char>= read_sequence().chars().collect();
        b.iter(|| SmithWaterman::new_thread(&g,&r, 32));
    }
    #[bench]
    fn bench_maxtrix_thread_many(b: &mut Bencher) {
        let g= genome_sequence_many(1_000).chars().collect();
        let r= read_sequence_many(1_000).chars().collect();
        b.iter(|| SmithWaterman::new_thread(&g,&r, 32));
    }
    #[bench]
    fn bench_maxtrix_thread_many_many(b: &mut Bencher) {
        let g= genome_sequence_many(10_000).chars().collect();
        let r= read_sequence_many(1_000).chars().collect();
        b.iter(|| SmithWaterman::new_thread(&g,&r, 32));
    }
}

