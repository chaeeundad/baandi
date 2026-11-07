//! Test-only source injected into an isolated measurement copy by the engine harness.
//! Run: cargo test --release -p dft-engine --test generalized_gram_stage0 -- --ignored --nocapture
//! No production module includes this file. The source string uses the already adopted
//! faer dependency, added only to the measurement copy's manifest.
pub const PROBE: &str = r###"
use crate::{BandTolerances, LinearOperator, RealLinearOperator, Preconditioner};
use dft_numerics::{linalg::{ComplexMatrix, hermitian_eigh}, real::{RealMatrix, symmetric_eigh}};
use faer::{Mat, Side};
use num_complex::Complex64 as C;
use std::sync::atomic::{AtomicUsize, Ordering};
static IDS: AtomicUsize = AtomicUsize::new(0);
fn cm(a: &Mat<C>) -> ComplexMatrix {
    ComplexMatrix::new(a.nrows(), a.ncols(), (0..a.ncols()).flat_map(|j| (0..a.nrows()).map(move |i| a[(i,j)])).collect()).unwrap()
}
fn rm(a: &Mat<C>) -> RealMatrix {
    RealMatrix::new(a.nrows(), a.ncols(), (0..a.ncols()).flat_map(|j| (0..a.nrows()).map(move |i| { assert_eq!(a[(i,j)].im,0.); a[(i,j)].re })).collect()).unwrap()
}
fn from_c(a: &ComplexMatrix) -> Mat<C> { Mat::from_fn(a.rows(), a.cols(), |i,j| a[(i,j)]) }
fn from_r(a: &RealMatrix) -> Mat<C> { Mat::from_fn(a.rows(), a.cols(), |i,j| C::new(a[(i,j)],0.)) }
fn eig(a: &Mat<C>, real: bool) -> (Vec<f64>, Mat<C>) {
    if real { let (s,u)=symmetric_eigh(&rm(a)).unwrap(); (s,from_r(&u)) }
    else { let (s,u)=hermitian_eigh(&cm(a)).unwrap(); (s,from_c(&u)) }
}
fn sym(a: &mut Mat<C>) {
    for j in 0..a.ncols() { a[(j,j)].im=0.; for i in 0..j { a[(j,i)]=a[(i,j)].conj(); } }
}
// faer's maintained kernels, explicitly sequential. Gamma products use f64.
fn mul(a: &Mat<C>, b: &Mat<C>, adj: bool, real: bool) -> Mat<C> {
    let m=if adj {a.ncols()} else {a.nrows()};
    if real {
        let ar=Mat::from_fn(a.nrows(),a.ncols(),|i,j|a[(i,j)].re);
        let br=Mat::from_fn(b.nrows(),b.ncols(),|i,j|b[(i,j)].re);
        let mut out=Mat::<f64>::zeros(m,b.ncols());
        faer::linalg::matmul::matmul(out.as_mut(),faer::Accum::Replace,if adj {ar.transpose()} else {ar.as_ref()},br.as_ref(),1.,faer::Par::Seq);
        Mat::from_fn(m,b.ncols(),|i,j|C::new(out[(i,j)],0.))
    } else {
        let mut out=Mat::<C>::zeros(m,b.ncols());
        if adj { faer::linalg::matmul::matmul(out.as_mut(),faer::Accum::Replace,a.adjoint(),b.as_ref(),C::new(1.,0.),faer::Par::Seq); }
        else { faer::linalg::matmul::matmul(out.as_mut(),faer::Accum::Replace,a.as_ref(),b.as_ref(),C::new(1.,0.),faer::Par::Seq); }
        out
    }
}
fn cholesky_whitener(s: &Mat<C>, real: bool) -> Option<Mat<C>> {
    let m=s.nrows();
    if real {
        let a=Mat::from_fn(m,m,|i,j|s[(i,j)].re);
        let factor=a.llt(Side::Lower).ok()?;
        let mut inv=Mat::<f64>::identity(m,m);
        factor.L().solve_lower_triangular_in_place(inv.as_mut());
        Some(Mat::from_fn(m,m,|i,j|C::new(inv[(j,i)],0.)))
    } else {
        let factor=s.llt(Side::Lower).ok()?;
        let mut inv=Mat::<C>::identity(m,m);
        factor.L().solve_lower_triangular_in_place(inv.as_mut());
        Some(Mat::from_fn(m,m,|i,j|inv[(j,i)].conj()))
    }
}
fn oracle() {
    for real in [true,false] {
        let v=Mat::from_fn(3,3,|i,j|C::new(if i==j{1.}else{0.15*(i+j+1) as f64},if real{0.}else{0.1*(i as f64-j as f64)}));
        let hv=Mat::from_fn(3,3,|i,j|v[(i,j)]*(i+1) as f64);
        let mut s=mul(&v,&v,true,real);sym(&mut s);
        let h=mul(&v,&hv,true,real);
        let w=cholesky_whitener(&s,real).unwrap();
        let hw=mul(&h,&w,false,real);let mut t=mul(&w,&hw,true,real);sym(&mut t);
        let (e,z)=eig(&t,real);let c=mul(&w,&z,false,real);
        let x=mul(&v,&c,false,real);let hx=mul(&hv,&c,false,real);
        for j in 0..3 {assert!((e[j]-(j+1) as f64).abs()<1e-12);for i in 0..3{assert!((hx[(i,j)]-x[(i,j)]*e[j]).norm()<1e-12);}}
        let singular=Mat::from_fn(2,2,|_,_|C::new(1.,0.));
        assert!(cholesky_whitener(&singular,real).is_none());
        let (e,_)=eig(&singular,real);assert!(e[0].abs()<1e-12);assert!((e[1]-2.).abs()<1e-12);
    }
    println!("P37A_ORACLE,real_and_complex,analytic_spectrum_and_residual_lt_1e-12_Ha,singular_cholesky_rejected");
}
fn cols(a:&Mat<C>, ids:&[usize])->Mat<C>{Mat::from_fn(a.nrows(),ids.len(),|i,j|a[(i,ids[j])])}
fn append(a:&Mat<C>, b:&Mat<C>)->Mat<C>{Mat::from_fn(a.nrows(),a.ncols()+b.ncols(),|i,j|if j<a.ncols(){a[(i,j)]}else{b[(i,j-a.ncols())]})}
#[allow(clippy::too_many_arguments)]
fn run(initial: Mat<C>, bands:usize, strict:usize, tol:f64, buffer:f64, ritz:bool, max:usize, real:bool,
       apply:impl Fn(&Mat<C>)->Mat<C>, pre:impl Fn(&[f64],&Mat<C>)->Mat<C>) -> usize {
    let id=IDS.fetch_add(1,Ordering::Relaxed);
    if id==0 {oracle();}
    for mult in [2,3,6] {
        let cap=(mult*bands).min(initial.nrows());
        let mut v=initial.clone(); let mut hv=apply(&v);
        let mut old=0; let mut h=Mat::<C>::zeros(0,0); let mut s=h.clone();
        let mut previous=vec![f64::NAN;bands];
        let mut finished=false;
        for step in 1..=max {
            let n=v.nrows(); let m=v.ncols(); let added=m-old;
            let new:Vec<_>=(old..m).collect();
            let sh=mul(&v,&cols(&v,&new),true,real);
            let hh=mul(&v,&cols(&hv,&new),true,real);
            s=Mat::from_fn(m,m,|i,j|if i<old&&j<old{s[(i,j)]}else if j>=old{sh[(i,j-old)]}else{sh[(j,i-old)].conj()});
            h=Mat::from_fn(m,m,|i,j|if i<old&&j<old{h[(i,j)]}else if j>=old{hh[(i,j-old)]}else{hh[(j,i-old)].conj()});
            sym(&mut s);sym(&mut h);
            // Independent diagnostic Gram rebuild, excluded from candidate work.
            // It measures the actual floating-point basis, including contraction drift.
            let mut actual_s=mul(&v,&v,true,real);sym(&mut actual_s);
            let gram_drift=(0..m).flat_map(|j|(0..m).map(move|i|(i,j))).map(|(i,j)|(actual_s[(i,j)]-s[(i,j)]).norm()).fold(0_f64,f64::max);
            s=actual_s;
            let (sv,su)=eig(&s,real); let largest=sv[m-1]; let rel=sv[0]/largest;
            let rank=sv.iter().filter(|&&x|x>1e-10*largest).count();
            let chol_w=cholesky_whitener(&s,real);
            let chol=chol_w.is_some();
            let cond=if sv[0]>0.{largest/sv[0]}else{f64::INFINITY};
            // Below 1e-8 (or Cholesky failure), continuation uses canonical whitening with a fixed 1e-12
            // floor, retaining the pre-fallback spectrum in every output row.
            // This is diagnostic continuation, never an SCF solver result.
            let keep:Vec<_>=(0..m).filter(|&j|sv[j]>1e-12*largest).collect();
            let rect=2_u64*n as u64*m as u64*added as u64;
            if keep.len()<bands {
                println!("P37A_GRAM,{id},{real},{ritz},{mult},{step},{n},{m},{rank},{:.17e},{rel:.17e},{cond:.17e},{chol},{rect},0,0,0,rank_exhausted",sv[0]);
                break;
            }
            let w=if rel>=1e-8 && chol {chol_w.unwrap()}else{Mat::from_fn(m,keep.len(),|i,j|su[(i,keep[j])]/sv[keep[j]].sqrt())};
            let hw=mul(&h,&w,false,real); let mut reduced=mul(&w,&hw,true,real);sym(&mut reduced);
            let (values,z)=eig(&reduced,real);let ids:Vec<_>=(0..bands).collect();
            let c=mul(&w,&cols(&z,&ids),false,real);
            let threshold=|j:usize|if j<strict-usize::from(ritz){tol}else{buffer};
            // Strict must inspect every residual; do not claim unconverged-only
            // vector savings for its convergence decision.
            let selected:Vec<_>=if ritz {(0..bands).filter(|&j|!previous[j].is_finite()||(values[j]-previous[j]).abs()>threshold(j)).collect()}else{ids.clone()};
            let cc=cols(&c,&selected);let x=mul(&v,&cc,false,real);let hx=mul(&hv,&cc,false,real);
            let r=Mat::from_fn(n,selected.len(),|i,j|hx[(i,j)]-x[(i,j)]*values[selected[j]]);
            let active:Vec<_>=(0..selected.len()).filter(|&j|ritz||((0..n).map(|i|r[(i,j)].norm_sqr()).sum::<f64>()).sqrt()>threshold(selected[j])).collect();
            let rot=2_u64*n as u64*m as u64*selected.len() as u64;
            let restart=m+active.len()>cap;
            let contraction=if restart||active.is_empty(){2_u64*n as u64*m as u64*bands as u64}else{0};
            // Explicit model, not a faer instruction count: two full EVDs and
            // whitening products actually used by this diagnostic.
            let k=keep.len() as u64;let mm=m as u64;
            let small=2*mm.pow(3)+2*k.pow(3)+mm*mm*k+mm*k*k+mm*k*bands as u64+if restart {2*(mm*mm*bands as u64+mm*(bands as u64).pow(2))}else{0};
            println!("P37A_GRAM,{id},{real},{ritz},{mult},{step},{n},{m},{rank},{:.17e},{rel:.17e},{cond:.17e},{chol},{rect},{rot},{contraction},{small},{},{gram_drift:.17e}",sv[0],if active.is_empty(){"converged"}else{"continue"});
            if active.is_empty(){finished=true;break;}
            let shifts:Vec<_>=active.iter().map(|&j|values[selected[j]]).collect();
            let mut correction=pre(&shifts,&cols(&r,&active));
            for j in 0..correction.ncols(){let norm=(0..n).map(|i|correction[(i,j)].norm_sqr()).sum::<f64>().sqrt(); assert!(norm.is_finite()&&norm>0.);for i in 0..n{correction[(i,j)]/=norm;}}
            if restart {
                let hc=mul(&h,&c,false,real);h=mul(&c,&hc,true,real);
                let sc=mul(&s,&c,false,real);s=mul(&c,&sc,true,real);
                v=mul(&v,&c,false,real);hv=mul(&hv,&c,false,real);
                old=bands;
            } else {old=m;}
            let hc=apply(&correction);v=append(&v,&correction);hv=append(&hv,&hc);
            previous.copy_from_slice(&values[..bands]);
        }
        println!("P37A_END,{id},{mult},{finished}");
    }
    id
}
pub fn real(op:&dyn RealLinearOperator, initial:&RealMatrix, diagonal:&[f64], bands:usize, policy:BandTolerances, tol:f64, ritz:bool, max:usize)->usize {
    run(from_r(initial),bands,policy.strict_bands,tol,policy.buffer_tolerance,ritz,max,true,
    |x|{let x=rm(x);let mut y=RealMatrix::zeros(x.rows(),x.cols()).unwrap();op.apply_block(&x,&mut y).unwrap();from_r(&y)},
    |shifts,r|Mat::from_fn(r.nrows(),r.ncols(),|i,j|{let d=diagonal[i]-shifts[j];-r[(i,j)]/if d<0.{-d.abs().max(0.1)}else{d.abs().max(0.1)}}))
}
pub fn complex(op:&dyn LinearOperator, initial:&ComplexMatrix, p:Option<&dyn Preconditioner>, bands:usize, policy:BandTolerances, tol:f64, ritz:bool, max:usize)->usize {
    run(from_c(initial),bands,policy.strict_bands,tol,policy.buffer_tolerance,ritz,max,false,
    |x|{let x=cm(x);let mut y=ComplexMatrix::zeros(x.rows(),x.cols()).unwrap();op.apply_block(&x,&mut y).unwrap();from_c(&y)},
    |shifts,r|{if let Some(p)=p {let r=cm(r);let mut y=ComplexMatrix::zeros(r.rows(),r.cols()).unwrap();p.apply(shifts,&r,&mut y).unwrap();from_c(&y)}else{r.clone()}})
}
"###;
