#![feature(test)]
extern crate test;
extern crate smith_waterman;
mod tests{
    use smith_waterman::*;
    #[test]
    fn it_sets_the_matrix_fn_no_match(){
        let g = "ab".chars().collect();
        let r = "cd".chars().collect();
        let mut smitty = SmithWaterman::new(&g,&r);
        smitty.set_matrix_fn();
        let alignment = smitty.align();
        assert_eq!(("".to_string(),"".to_string()), alignment);
    }
    #[test]
    fn it_sets_the_matrix_loop_no_match(){
        let g = "ab".chars().collect();
        let r = "cd".chars().collect();
        let mut smitty = SmithWaterman::new(&g,&r);
        smitty.set_matrix_loops();
        let alignment = smitty.align();
        assert_eq!(("".to_string(),"".to_string()), alignment);
    }

    #[test]
    fn it_sets_the_matrix_fn_one_match(){
        let g = "ab".chars().collect();
        let r = "cb".chars().collect();
        let mut smitty = SmithWaterman::new(&g,&r);
        smitty.set_matrix_fn();
        assert_eq!(2, smitty.matrix[(2,2)]);
    }
    #[test]
    fn it_sets_the_matrix_loop_one_match(){
        let g = "ab".chars().collect();
        let r = "cb".chars().collect();
        let mut smitty = SmithWaterman::new(&g,&r);
        smitty.set_matrix_loops();
        assert_eq!(2, smitty.matrix[(2,2)]);
    }
    #[test]
    fn it_sets_aligns_wiki_example(){
        let g = "ACACACTA".chars().collect();
        let r = "AGCACACA".chars().collect();
        let mut smitty = SmithWaterman::new(&g,&r);
        let alignment = smitty.align();
        assert_eq!(("A-CACACTA".to_string(), "AGCACAC-A".to_string()), alignment);
    }
}

