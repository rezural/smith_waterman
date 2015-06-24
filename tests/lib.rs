#![feature(test)]

extern crate test;
extern crate smith_waterman;
mod tests{
    use smith_waterman::*;
    use test::Bencher;

    #[test]
    fn it_sets_the_matrix_fn_no_match(){
        let mut smitty = SmithWaterman::new("ab".to_string(), "cd".to_string());
        smitty.set_matrix_fn();
        let alignment = smitty.align();
        assert_eq!(("".to_string(),"".to_string()), alignment);
        assert!(smitty.matrix.is_zero());
    }
    #[test]
    fn it_sets_the_matrix_loop_no_match(){
        let mut smitty = SmithWaterman::new("ab".to_string(), "cd".to_string());
        smitty.set_matrix_loops();
        let alignment = smitty.align();
        assert_eq!(("".to_string(),"".to_string()), alignment);
        assert!(smitty.matrix.is_zero());
    }

    #[test]
    fn it_sets_the_matrix_fn_one_match(){
        let mut smitty = SmithWaterman::new("ab".to_string(), "cb".to_string());
        smitty.set_matrix_fn();
        assert_eq!(2, smitty.matrix[(2,2)]);
    }
    #[test]
    fn it_sets_the_matrix_loop_one_match(){
        let mut smitty = SmithWaterman::new("ab".to_string(), "cb".to_string());
        smitty.set_matrix_loops();
        assert_eq!(2, smitty.matrix[(2,2)]);
    }
    #[test]
    fn it_sets_aligns_wiki_example(){
        let mut smitty = SmithWaterman::new( "ACACACTA".to_string(),"AGCACACA".to_string());
        let alignment = smitty.align();
        assert_eq!(("A-CACACTA".to_string(), "AGCACAC-A".to_string()), alignment);
    }

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
        let mut smitty = SmithWaterman::new(genome_sequence(), read_sequence());
        b.iter(|| smitty.set_matrix_loops());
    }

    #[bench]
    fn bench_maxtrix_loops_many(b: &mut Bencher) {
        return;
        let mut smitty = SmithWaterman::new(genome_sequence_many(1_000), read_sequence_many(1_000));
        b.iter(|| smitty.set_matrix_loops());
    }
    #[bench]
    fn bench_maxtrix_loops_many_many(b: &mut Bencher) {
        return;
        let mut smitty = SmithWaterman::new(genome_sequence_many(10_000), read_sequence_many(10_000));
        b.iter(|| smitty.set_matrix_loops());
    }
    #[bench]
    fn bench_maxtrix_fn(b: &mut Bencher) {
        let mut smitty = SmithWaterman::new(genome_sequence(), read_sequence());
        b.iter(|| smitty.set_matrix_fn());
    }
    #[bench]
    fn bench_maxtrix_fn_many(b: &mut Bencher) {
        return;
        let mut smitty = SmithWaterman::new(genome_sequence_many(1_000), read_sequence_many(1_000));
        b.iter(|| smitty.set_matrix_fn());
    }
    #[bench]
    fn bench_maxtrix_fn_many_many(b: &mut Bencher) {
        return;
        let mut smitty = SmithWaterman::new(genome_sequence_many(10_000), read_sequence_many(10_000));
        b.iter(|| smitty.set_matrix_fn());
    }
}

