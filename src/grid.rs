use std::sync::Arc;
use polars::prelude::{DataFrame, CsvReader, SerReader, PolarsError, DataType, Field, Schema};


#[derive(Debug)]
pub struct GRanges {
    pub seqnames: Vec<String>,
    pub start: Vec<i64>,
    pub end: Vec<i64>,
}


#[derive(Debug)]
pub struct Block {
    pub gro: Vec<GRanges>,
    pub cell: Vec<String>,
    pub sv_state: Vec<String>
}

#[derive(Debug)]
pub struct Grid {
    pub gro: Vec<GRanges>,
    pub id: Vec<i64>
}

pub fn read_sv_data() -> Result<DataFrame, PolarsError> {

    let schema = Schema::from_iter(vec![
        Field::new("llr_to_ref", DataType::Float64)
    ]);

    let df = CsvReader::from_path("data/tall3_proc.tsv.gz")
        .expect("cannot open file")
        .has_header(true)
        .with_separator(b'\t')
        .with_dtypes(Some(Arc::new(schema)))
        .finish();

    return df;
}
