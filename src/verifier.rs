

use crate::prover;

pub use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    polynomial::univariate::SparsePolynomial as UniSparsePolynomial,
    DenseMVPolynomial, Polynomial};

pub use ark_bls12_381::Fq;



pub struct Verifier {
    pub g: prover::MultiPoly, 
    pub v: Vec<i32>
}

impl Verifier {

    pub fn new(g: &prover::MultiPoly, v: Vec<i32>) -> Self { 
        Verifier {
            g: g.clone(),
            v: v.clone()
	    }
    }

    pub fn verify_sum_check_protocol (&mut self) -> bool { 
        let mut p = prover::Prover::new(&self.g);
        let mut c1 = p.compute_h(); 

        println!("H: {}", cast_fq_to_i32(c1));
        let mut r_vec: Vec<Fq> = vec![]; 
        let mut si: prover::UniPoly;
        let mut ok = true;
        for i in 0..self.g.num_vars() {
            si = p.compute_sum(r_vec.clone(), i);
            print!("s{} = ", i+1);
            print_unipol(si.clone());

            if verifier_condition(si.clone(), c1) {
                let r = get_r(i, self.v.clone());
                r_vec.insert(r_vec.len(),r);
                c1 = evaluate_gi(si.clone(), r);
                
                continue;
            }
            ok = false;
            break;
        }

        if ok == true {
            let last = Polynomial::evaluate(&self.g.clone(), &r_vec);
            if last == c1 {
                println!("Sum Check verification passed!");
                return true;
            }
        }
        
        println!("Sum Check verification NOT passed!");
        return false;
    }
}

fn evaluate_gi (gi: prover::UniPoly, value: Fq) -> Fq {
    let s = Polynomial::evaluate(&gi.clone(), &value.clone());
    println!("s[{}] = {}", value, s);
    s
}

fn verifier_condition (gi: prover::UniPoly, si: Fq ) -> bool {
    evaluate_gi(gi.clone(), Fq::from(0))+ evaluate_gi(gi.clone(), Fq::from(1)) == si
}


fn get_r(i: usize, random_vec: Vec<i32>) -> Fq{
    //let mut rng = rand::thread_rng();
	// r :Fq = rng.gen();

	Fq::from(random_vec[i])
}

fn cast_fq_to_i32 (f: Fq) -> i32 {
    f.to_string().parse::<i32>().unwrap()
}

fn print_unipol(si: prover::UniPoly) {
    print!("{:?} \n", (for i in si.clone().to_vec() {
        print!("{}*x^{}+", cast_fq_to_i32(i.1), i.0);
    }));
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sumcheck_test() {
        let g: prover::MultiPoly = SparsePolynomial::from_coefficients_vec(
                4,
                vec![
                    (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
                    (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
                    (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),  
                ],
            );
        let r_random_vector: Vec<i32> = vec![2,3,6];

        let mut v = Verifier::new(&g, r_random_vector);

        assert!(v.verify_sum_check_protocol()); 
    }
}