#![feature(test)]
extern crate test;
extern crate smith_waterman;
mod tests{
    use smith_waterman::*;
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
}

