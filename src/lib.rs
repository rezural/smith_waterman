#![feature(str_char)]
extern crate nalgebra;

use nalgebra::DMat;
use std::fmt::{Debug, Formatter, Result};

pub struct SmithWaterman{
    genome_sequence: String,
    read_sequence: String,
    matrix:  DMat<isize>,
    matched: isize,
    missed: isize,
}
pub enum GraphMovements {Blank, Left, Top, Diagonal}

impl SmithWaterman{

    pub fn new(genome_sequence: String, read_sequence: String) -> SmithWaterman {
        SmithWaterman{matrix: nalgebra::DMat::new_zeros(read_sequence.len()+1, genome_sequence.len()+1),
        genome_sequence: genome_sequence, read_sequence: read_sequence, matched: 2, missed: -1}
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
        for row in (1..self.read_sequence.len()+1){
            for col in (1..self.genome_sequence.len()+1){
                let left = self.penalty(self.matrix[(row, col-1)], self.missed);
                let top = self.penalty(self.matrix[(row-1, col)], self.missed);
                let diagonal = self.matrix[(row-1, col-1)];
                let diagonal_match = if self.read_sequence.char_at(row-1) == self.genome_sequence.char_at(col-1){
                    self.penalty(diagonal, self.matched)
                }else{
                    self.penalty(diagonal, self.missed)
                };
                let n = std::cmp::max(left, std::cmp::max(top, diagonal_match));
                if n >= max{
                    max = n;
                    max_point = (row,col);
                }
                self.matrix[(row, col)] = n;
            }
        }
        let mut last_movement = GraphMovements::Blank;
        let mut genome_sequence_alignment = String::new();
        let mut read_sequence_alignment = String::new();
        while max > 0 {
            let (row, col) = max_point;
            let one = self.genome_sequence.char_at(col-1);
            let two = self.read_sequence.char_at(row-1);
            match last_movement{
                GraphMovements::Blank  => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push(two);
                },
                GraphMovements::Diagonal => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push(two);
                },
                GraphMovements::Top => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push('-');
                },
                GraphMovements::Left => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push('-');
                },
            }
            let top = self.matrix[(row-1, col)];
            let left = self.matrix[(row, col-1)];
            let diagonal = self.matrix[(row-1, col-1)];
            max  = std::cmp::max(left, std::cmp::max(top, diagonal));
            max_point = if diagonal == max{
                last_movement = GraphMovements::Diagonal;
                (row-1, col-1)
            } else if left == max{
                last_movement = GraphMovements::Left;
                (row, col-1)
            } else {
                last_movement = GraphMovements::Top;
                (row-1, col)
            };
        };
        let x1: String = genome_sequence_alignment.chars().rev().collect();
        let x2: String = read_sequence_alignment.chars().rev().collect();
        println!("{:?}",x1);
        println!("{:?}",x2);
    }
}

impl Debug for SmithWaterman {
    fn fmt(&self, form:&mut Formatter) -> Result {
        //nrows already has an extra over the sequence counts for the row of zeros
        for row in 0..self.matrix.nrows()+1 {
            for col in 0..self.matrix.ncols()+1 {
                let _ = if col==0 && row>1{
                    write!(form, "{:>5}", self.read_sequence.char_at(row-2).to_string())
                } else if row==0 && col>1{
                    write!(form, "{:>5}", self.genome_sequence.char_at(col-2).to_string())
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
