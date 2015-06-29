extern crate nalgebra;
use nalgebra::DMat;
use std::fmt::{Debug, Formatter, Result};

use std::sync::Arc;
use std::sync::mpsc::{channel, Sender, Receiver};
use std::thread;
use std::sync::mpsc;
use std::collections::HashSet;


/// The `SmithWaterman` struct
///
/// genome_sequence: String
///
/// read_sequence: String
///
/// matrix:  DMat<isize>
pub struct SmithWaterman {
    /// Genome
    pub genome_sequence: Vec<char>,
    /// Compared Genome
    pub  read_sequence:  Vec<char>,
    /// Matrix used to store data
    pub matrix:  DMat<isize>,
    pub matched: isize,
    pub missed: isize,
}
pub enum GraphMovements {Blank, Left, Top, Diagonal}

impl SmithWaterman{
    /// Constructs a new `SmithWaterman`.
    ///
    /// # Examples
    ///
    /// ```
    /// use smith_waterman;
    ///
    /// let mut smitty = smith_waterman::SmithWaterman::new("ab".to_string(), "cb".to_string());
    /// ```
    pub fn new(genome_sequence: String, read_sequence: String) -> SmithWaterman {
        let matrix = DMat::new_zeros(0,0);
        SmithWaterman{matrix: matrix, genome_sequence: genome_sequence.chars().collect(),
        read_sequence: read_sequence.chars().collect(), matched: 2, missed: -1}
    }

    /// Constructs a new `SmithWaterman`.
    ///
    /// # Examples
    ///
    /// ```
    /// use smith_waterman;
    ///
    /// let mut smitty = smith_waterman::SmithWaterman::new_chars("ab".chars().collect(), "cb".chars().collect());
    /// ```
    pub fn new_chars(genome_sequence: Vec<char>, read_sequence: Vec<char>) -> SmithWaterman {
        let matrix = DMat::new_zeros(0,0);
        SmithWaterman{matrix: matrix, genome_sequence: genome_sequence,
        read_sequence: read_sequence, matched: 2, missed: -1}
    }

    /// Constructs and fills a new `SmithWaterman`.
    ///
    /// # Examples
    ///
    /// ```
    /// use smith_waterman;
    ///
    /// let mut smitty = smith_waterman::SmithWaterman::new_threads("ab".chars().collect(), "cb".chars().collect());
    /// ```
    pub fn new_threads(genome_sequence: Vec<char>, read_sequence: Vec<char>) -> SmithWaterman {
        let matrix = nalgebra::DMat::new_zeros(read_sequence.len()+1, genome_sequence.len()+1);
        let shared_matrix = Arc::new(matrix);
        let shared_genome = Arc::new(genome_sequence);
        let shared_read = Arc::new(read_sequence);
        let mut max_point = (0,0);
        let mut max = 0;
        let thread_counter = 0;
        let seen_locations: HashSet<(usize, usize)> = HashSet::new();
        let (fromParentSender, fromParentReceiver) = channel();
        let (fromChildSender, fromChildReceiver) = channel();
        for id in 0..2 {
            let thread_sender = fromChildSender.clone();
            let thread_matrix = shared_matrix.clone();
            let thread_genome = shared_genome.clone();
            let thread_read = shared_read.clone();
            let guard = thread::spawn(|| {
                for thread_data in fromParentReceiver{
                    let (row, col) = (0,9);//thread_data;
                    let matched = thread_genome.get(row-1).unwrap() ==
                thread_read.get(col-1).unwrap();
                    let n = SmithWaterman::threaded_calculated_movement(&thread_matrix, matched,  row, col, 0-1, 2);
                    thread_sender.send(1);
                }
            });
        }

        fromParentSender.send(1);

        loop{
            let x = fromChildReceiver.recv().unwrap();
            let (row, col, value) = (2,3,4);
            if value >= max{
                max = value;
                max_point = (row, col);
            }
            matrix[(row, col)] = value;
            break;
        }


        SmithWaterman{matrix: matrix, genome_sequence: genome_sequence,
        read_sequence: read_sequence, matched: 2, missed: 0-1}
    }

    fn threaded_calculated_movement(matrix: &DMat<isize>, matches: bool, row: usize, col: usize, missed: isize, matched: isize) ->
        isize{
        let left = if col>=1 {SmithWaterman::penalty_thread(matrix[(row, col-1)], missed)} else {0};
        let top = if row>=1 {SmithWaterman::penalty_thread(matrix[(row-1, col)], missed)} else {0};
        let diagonal_value = if row>=1 && col>=1 {matrix[(row-1, col-1)]} else {0};
        let diagonal_match = if row == 0 || col == 0{
            0
        }else if matches{
            SmithWaterman::penalty_thread(diagonal_value, matched)
        }else{
            SmithWaterman::penalty_thread(diagonal_value, missed)
        };
        std::cmp::max(left, std::cmp::max(top, diagonal_match))
    }

    fn penalty_thread(value: isize, penalty_value: isize) -> isize{
        match value.checked_add(penalty_value){
            Some(i) =>{
                if i<0 { 0 }else { i }
            },
            _ => {0}
        }
    }

    fn penalty(&self, value: isize, penalty_value: isize) -> isize{
        match value.checked_add(penalty_value){
            Some(i) =>{
                if i<0 { 0 }else { i }
            },
            _ => {0}
        }
    }

    /// Runs the matrix to obtain the optimum local alignment.
    /// Uses the set_matrix_loops function.
    /// returns a tuple of strings. The fist is the genome sequence the second is the read
    /// sequence.
    /// Faster in smaller sets: https://gist.github.com/sbeckeriv/d9d2c03b19178a888c32
    /// # Examples
    ///
    /// ```
    /// let mut smitty = smith_waterman::SmithWaterman::new("ab".to_string(), "cb".to_string());
    /// let alignment = smitty.align();
    /// ```
    pub fn align(&mut self) -> (String, String){
        let max_point = self.set_matrix_loops();
        return self.local_alignment(max_point);
    }

    /// Runs the matrix to obtain the optimum local alignment.
    /// Uses the set_matrix_fn function.
    /// returns a tuple of strings. The fist is the genome sequence the second is the read
    /// sequence.
    ///
    /// In benchmarks this function can be 2x faster. Tested with a 70k and 7k string. <Caveat
    /// Emptor>
    /// https://gist.github.com/sbeckeriv/d9d2c03b19178a888c32
    /// # Examples
    ///
    /// ```
    /// let mut smitty = smith_waterman::SmithWaterman::new("ab".to_string(), "cb".to_string());
    /// let alignment = smitty.align_fn();
    /// ```
    pub fn align_fn(&mut self) -> (String, String){
        let max_point = self.set_matrix_fn();
        return self.local_alignment(max_point);
    }

    /// Fills the matrix with values.
    /// This uses a single loop and calculates the correct col and row on the fly.
    /// returns a tuple of usize. The fist is the row and the second is the column of the largest
    /// value in the matrix.
    /// This value can also be calculated out from `smitty.matrix`
    /// In benchmarks this function can be 2x faster. Tested with a 70k and 7k string. <Caveat
    /// Emptor>
    /// https://gist.github.com/sbeckeriv/d9d2c03b19178a888c32
    ///
    /// # Examples
    /// ```
    /// let mut smitty = smith_waterman::SmithWaterman::new("ab".to_string(), "cb".to_string());
    /// let max_point = smitty.set_matrix_fn();
    /// ```
    pub fn set_matrix_fn(&mut self) -> (usize, usize) {
        let mut max_point = (0,0);
        let mut max = 0;
        let rows = self.read_sequence.len()+1;
        let cols = self.genome_sequence.len()+1;
        self.matrix = nalgebra::DMat::new_zeros(self.read_sequence.len()+1, self.genome_sequence.len()+1);
        for i in (0..rows * cols){
            let col = i / rows;
            let row = i - col * rows;
            let n = self.calculated_movement(row, col);
            if n >= max{
                max = n;
                max_point = (row,col);
            }
            self.matrix[(row, col)]=n;
        }
        return max_point;
    }

    fn calculated_movement(&self, row: usize, col: usize) -> isize{
        let left = if col>=1 {self.penalty(self.matrix[(row, col-1)], self.missed)} else {0};
        let top = if row>=1 {self.penalty(self.matrix[(row-1, col)], self.missed)} else {0};
        let diagonal_value = if row>=1 && col>=1 {self.matrix[(row-1, col-1)]} else {0};
        let diagonal_match = if row == 0 || col == 0{
            0
        }else if self.read_sequence.get(row-1).unwrap() == self.genome_sequence.get(col-1).unwrap(){
            self.penalty(diagonal_value, self.matched)
        }else{
            self.penalty(diagonal_value, self.missed)
        };
        std::cmp::max(left, std::cmp::max(top, diagonal_match))
    }

    /// Fills the matrix with values.
    /// This uses a naïve double loop to fill the matrix.
    ///
    /// returns a tuple of usize. The fist is the row and the second is the column of the largest
    /// value in the matrix.
    /// This value can also be calculated out from `smitty.matrix`
    ///
    /// Faster in smaller sets: https://gist.github.com/sbeckeriv/d9d2c03b19178a888c32
    /// # Examples
    ///
    /// ```
    /// let mut smitty = smith_waterman::SmithWaterman::new("ab".to_string(), "cb".to_string());
    /// let max_point = smitty.set_matrix_loops();
    /// ```
    pub fn set_matrix_loops(&mut self) -> (usize, usize) {
        self.matrix = nalgebra::DMat::new_zeros(self.read_sequence.len()+1, self.genome_sequence.len()+1);
        let mut max_point = (0,0);
        let mut max = 0;
        for row in (1..self.read_sequence.len()+1){
            for col in (1..self.genome_sequence.len()+1){
                let n = self.calculated_movement(row, col);
                if n >= max{
                    max = n;
                    max_point = (row,col);
                }
                self.matrix[(row, col)] = n;
            }
        }
        return max_point;
    }

    fn local_alignment(&self, max_point_value: (usize, usize)) -> (String, String){
        let mut max_point = max_point_value;
        let mut max = self.matrix[max_point];
        let mut last_movement = GraphMovements::Blank;
        let mut genome_sequence_alignment: Vec<char> =  Vec::with_capacity(self.genome_sequence.len());
        let mut read_sequence_alignment: Vec<char> =  Vec::with_capacity(self.read_sequence.len());
        while max > 0 {
            let (row, col) = max_point;
            let one = self.genome_sequence.get(col-1).unwrap().clone();
            let two = self.read_sequence.get(row-1).unwrap().clone();
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

            match last_movement {
                GraphMovements::Blank  => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push(two);
                },
                GraphMovements::Diagonal => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push(two);
                },
                GraphMovements::Top => {
                    genome_sequence_alignment.push('-');
                    read_sequence_alignment.push(two);
                },
                GraphMovements::Left => {
                    genome_sequence_alignment.push(one);
                    read_sequence_alignment.push('-');
                },
            }
        };
        genome_sequence_alignment.reverse();
        let x1: String = genome_sequence_alignment.into_iter().collect();
        read_sequence_alignment.reverse();
        let x2: String = read_sequence_alignment.into_iter().collect();
        return (x1,x2)
    }
}

impl Debug for SmithWaterman {
    fn fmt(&self, form:&mut Formatter) -> Result {
        //nrows already has an extra over the sequence counts for the row of zeros
        for row in 0..self.matrix.nrows()+1 {
            for col in 0..self.matrix.ncols()+1 {
                let _ = if col==0 && row>1{
                    write!(form, "{:>5}", self.read_sequence.get(row-2).unwrap().to_string())
                } else if row==0 && col>1{
                    write!(form, "{:>5}", self.genome_sequence.get(col-2).unwrap().to_string())
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
    let alignment = smitty.align();
    println!("{:?}", smitty);
    println!("{:?}", alignment);

    let mut smitty =
        SmithWaterman::new_chars("ACACACTA".chars().collect(),"AGCACACA".chars().collect());
    let alignment = smitty.align_fn();
    println!("fn::{:?}", smitty);
    println!("fn::{:?}", alignment);

}
