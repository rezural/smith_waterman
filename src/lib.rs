#![feature(str_char)]
extern crate nalgebra;

use nalgebra::DMat;
use std::fmt::{Debug, Formatter, Result};

pub struct SmithWaterman{
    sequence1: String,
    sequence2: String,
    matrix:  DMat<isize>,
    matched: isize,
    missed: isize,
}

impl SmithWaterman{
    pub fn new(sequence1: String, sequence2: String) -> SmithWaterman {
        SmithWaterman{matrix: nalgebra::DMat::new_zeros(sequence1.len()+1, sequence2.len()+1),
            sequence1: sequence1, sequence2: sequence2, matched: 2, missed: -1}
    }

    fn penalty(&self, value: isize, penalty_value: isize) -> isize{
        match value.checked_add(penalty_value){
            Some(i) =>{
                if i<0 { 0 }else { i }
            },
            _ => {0}
        }
    }

    pub fn score(&mut self){
        let mut max_point = (0,0);
        let mut max = 0;
        for row in (1..self.sequence1.len()+1){
            for col in (1..self.sequence2.len()+1){
                let left = self.penalty(self.matrix[(row, col-1)], self.missed);
                let top = self.penalty(self.matrix[(row-1, col)], self.missed);
                let diagonal = self.matrix[(row-1, col-1)];
                let diagonal_match = if self.sequence2.char_at(col-1) == self.sequence1.char_at(row-1){
                    self.penalty(diagonal, self.matched)
                }else{
                    self.penalty(diagonal, self.missed)
                };
                let n = std::cmp::max(left, std::cmp::max(top, diagonal_match));
                if n>=max{
                    max =n;
                    max_point = (row,col);
                }
                self.matrix[(row, col)] = n;
            }
        }
    }
}
impl Debug for SmithWaterman {
    fn fmt(&self, form:&mut Formatter) -> Result {
        //nrows already has an extra over the sequence counts for the row of zeros
        for row in 0..self.matrix.nrows()+1 {
            for col in 0..self.matrix.ncols()+1 {
                let _ = if col==0 && row>1{
                    write!(form, "{:>5}", self.sequence1.char_at(row-2).to_string())
                } else if row==0 && col>1{
                    write!(form, "{:>5}", self.sequence2.char_at(col-2).to_string())
                } else if row>=1 && col>=1{
                    write!(form, "{:>5}", self.matrix[(row-1,col-1)])
                }else{
                    write!(form, "{:>5}", "-")
                };
            }
            let _ = write!(form, "\n");
        }
        write!(form, "\n")
    }
}

#[test]
fn its_debugging() {
    let mut smitty = SmithWaterman::new("atgcatgcatgc".to_string(), "atgggcatg".to_string());
    smitty.score();
    println!("{:?}", smitty);

    let mut smitty = SmithWaterman::new( "atgggcatg".to_string(),"atgcatgcatgc".to_string());
    smitty.score();
    println!("{:?}", smitty);
}
