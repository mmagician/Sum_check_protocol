use ark_bls12_381::Fq;

pub use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm},
    polynomial::univariate::SparsePolynomial as UniSparsePolynomial,
    DenseMVPolynomial, Polynomial};

pub type MultiPoly = SparsePolynomial<Fq, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<Fq>;

pub struct Prover {
    pub g: MultiPoly
}

impl Prover {

    pub fn new(g: &MultiPoly) -> Self { 
        Prover {
		g: g.clone()
	    }
    }

    pub fn compute_h(&mut self) -> Fq {
        let number_of_variants = self.g.num_vars;
        let number_of_combinations = 2u64.pow(number_of_variants as u32);
        let mut h = Fq::from(0);
        (0..number_of_combinations).for_each(|i| {
            let combination = get_fq_binary_representation(i, number_of_variants);
            h = h + Polynomial::evaluate(&self.g.clone(), &combination);
        });
        h
    }

    pub fn compute_sum(&mut self, r_vec: Vec<Fq>, index: usize) -> UniPoly {
        let unfixed_variants = self.g.num_vars - index - 1;
        
        let mut si: UniPoly = UniSparsePolynomial::default();
        let number_of_combinations = 2_u64.pow(unfixed_variants as u32);
    
        for i in 0..number_of_combinations {
            let mut gi = r_vec.clone();
            let mut combination_vector: Vec<Fq> = get_fq_binary_representation(i, unfixed_variants);
            gi.append(&mut combination_vector);
    
            //mark witk -1 the index position
            gi.insert(index, Fq::from(-1));
    
            si = si + compute_unipolinom(self.g.clone(), gi, index);
        }
        si
    }
 
    

}

fn compute_unipolinom(
    polynom: MultiPoly,
    values: Vec<Fq>,
    index: usize
) -> UniPoly {
    let mut coefficients_vector = vec![(0, Fq::from(0))];
    for terms in &polynom.terms {

        let mut zeroes = false;
        let mut variant = false;
        let mut variant_pow = 0;

        let mut current_terms_evaluation : Fq = Fq::from(terms.0);

        for term in terms.1.iter() {
            if term.0 == index {
                variant = true;
                variant_pow = term.1;
                continue;
            }

            if values[term.0] == Fq::from(0) {
                zeroes = true;
                break;
            }
            
            let v = cast_fq_to_i32(values[term.0]).pow(term.1 as u32);
            current_terms_evaluation = current_terms_evaluation * Fq::from(v);
        }

        if zeroes {
            continue; // if one vale of the multiplication terms is 0 => nothing to be added
        }
        if variant {
            //if x^variant_pow already exists in the coefficient vector
            if coefficients_vector.len() > variant_pow {
                //add current_terms_evaluation to existing coefficient
                coefficients_vector[variant_pow].1 = coefficients_vector[variant_pow].1 + current_terms_evaluation;
            }
            else {
                complete_v_with_0(&mut coefficients_vector, variant_pow);
                coefficients_vector.insert(variant_pow, (variant_pow, current_terms_evaluation));
            }
        }
        else {
            //add current_terms_evaluation to existing x^0 coeficients
            coefficients_vector[0].1 = coefficients_vector[0].1 + current_terms_evaluation;
        }   
    }

    UniPoly::from_coefficients_vec(coefficients_vector) 
}



fn cast_fq_to_i32 (f: Fq) -> i32 {
    f.to_string().parse::<i32>().unwrap()
}

fn complete_v_with_0 (v: &mut Vec<(usize, Fq)>, new_length: usize) {
    for i in v.len()..new_length {
        v.insert(i, (i, Fq::from(0)));
    }
}

fn get_fq_binary_representation(number: u64, pow: usize) -> Vec<Fq> {
    return (0..pow).map(|n| Fq::from((number >> n) & 1)).collect();
}