pub mod grid;
use grid::Grid;
use grid::GRanges;
use grid::Block;
use polars::lazy::dsl::col;
use polars::prelude::*;

fn main() {
    let grid = Grid {
        gro: vec![GRanges { seqnames: vec!["chr1".to_string()], start: vec![100], end: vec![150] }],
        id: vec![1]
    };

    let block = Block{
        gro: vec![GRanges { seqnames: vec!["chr1".to_string()], start: vec![120], end: vec![130] }],
        cell: vec!["cell1".to_string()],
        sv_state: vec!["amp".to_string()]
    };

    if block.gro[0].start > grid.gro[0].start {
        println!("Grid: {:?}\nBlock {:?}", grid, block);
    }

    let mut df = grid::read_sv_data()
        .expect("Could not read the df");

    df = df
        .lazy()
        .filter((col("bin_id").gt(168)).and(col("bin_id").lt(170)))
        .select([col("llr_to_ref"), col("chrom"), col("start_loc"), col("end_loc"), col("sv_call_name")])
        .collect()
        .expect("Cannot perform the operation");


    println!("{:?}", df);
}
