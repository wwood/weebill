use clap::Parser;
use sylph::cmdline::*;
use sylph::contain;
use sylph::inspect;
use sylph::merge;
use sylph::refdelta;
use sylph::sketch;
use sylph::twostage_db;
//use std::panic::set_hook;

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though.
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    //    set_hook(Box::new(|info| {
    //        if let Some(s) = info.payload().downcast_ref::<String>() {
    //            log::error!("{}", s);
    //        }
    //    }));
    let cli = Cli::parse();
    match cli.mode {
        Mode::Sketch(sketch_args) => sketch::sketch(sketch_args),
        Mode::Query(contain_args) => contain::contain(contain_args, false),
        Mode::Profile(contain_args) => contain::contain(contain_args, true),
        Mode::Inspect(inspect_args) => inspect::inspect(inspect_args),
        Mode::Merge(args) => merge::merge(args),
        Mode::RefBuild(args) => refdelta::run_ref_build(args),
        Mode::RefCompress(args) => refdelta::run_ref_compress(args),
        Mode::DbConvert(args) => twostage_db::run_db_convert(args),
    }
}
