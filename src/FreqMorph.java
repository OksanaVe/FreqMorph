package beast.evolution.substitutionmodel;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.StandardData;
import beast.evolution.tree.Node;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.substitutionmodel.Frequencies;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;


@Description ("Morph-Model as implemented by A.Gavryushkina & J.Heled with modifications to calculate empirical character state frequencies")
@Citation("Lewis, Paul O. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic biology 50.6(2001): 913 - 925.")

public class FM extends SubstitutionModel.Base {

	public Input<Integer> nrOfStatesInput = new Input<Integer>("stateNumber", "the number of character states");
	public Input<DataType> dataTypeInput = new Input<DataType>("datatype", "datatype, used to determine the number of states", Validate.XOR, nrOfStatesInput);
	public Input<Boolean> estimateInput = new Input<>("estimate", "whether to estimate the frequencies from data or assume uniform distribution over characters", true);
	public Input<TaxonSet> taxonSetInput = new Input<>("taxa", "An optional taxon-set used only to sort the sequences into the same order as they appear in the taxon-set.", new TaxonSet(), Validate.OPTIONAL);
    public Input<Alignment> dataInput = new Input<>("data", "Sequence data for which frequencies are calculated");

	boolean hasFreqs;
	private boolean updateFreqs;
	
	public FM () {
		frequenciesInput.setRule(Validate.OPTIONAL);
		try {
			
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("initAndValidate() call failed when constructing FM()");
		}
	}
	
	double totalSubRate;
	double[] frequencies; 
	int pattern;
	int counts;
	EigenDecomposition eigenDecomposition; 
	
	
	private void setFrequencies() {

		final Frequencies frequencies1 = frequenciesInput.get();
	    Frequencies frequencies2 = new Frequencies();
		if (frequencies1 != null) {
			if (frequencies1.getFreqs().length != nrOfStates) {
				throw new RuntimeException ("number of stationary frequencies does not match number of states.");
			}
			System.arraycopy(frequencies1.getFreqs(), 0, frequencies, 0, nrOfStates);
			totalSubRate = 1;
			for (int k = 0; k < nrOfStates; ++k) {
				totalSubRate -= frequencies[k]*frequencies[k];
			}
			hasFreqs = true; 
		} else {
		      frequencies = frequencies2.getFreqs();
		    	 
		      }
		Log.info.println("Starting frequencies: " + Arrays.toString(frequencies));
		 	hasFreqs = false;
	}
	  
	@Override
	public void initAndValidate() {
		if (nrOfStatesInput.get() != null) {
			nrOfStates = nrOfStatesInput.get();
		} else {
			nrOfStates = dataTypeInput.get().getStateCount();
			
		}
        frequencies = new double[nrOfStates];
        setFrequencies();     
      
      
    }
	
	
	
	
	 @Override
	    public double[] getFrequencies() {
	        return frequencies;
	    }

	    @Override
	    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
	        if( updateFreqs ) {
	           setFrequencies();
	        }

	        if( hasFreqs ) {
	            final double e1 = Math.exp(-(fStartTime - fEndTime) * fRate/totalSubRate);
	            final double e2 = 1 - e1;

	            for( int i = 0; i < nrOfStates; ++i ) {
	                final int r = i * nrOfStates;
	                for( int j = 0; j < nrOfStates; ++j ) {
	                    matrix[r + j] = frequencies[j] * e2;
	                }
	                matrix[r + i] += e1;
	            }
	        } else {
	            double fDelta = (nrOfStates / (nrOfStates - 1)) * (fStartTime - fEndTime);
	            double fPStay = (1.0 + (nrOfStates - 1) * Math.exp(-fDelta * fRate)) / nrOfStates;
	            double fPMove = (1.0 - Math.exp(-fDelta * fRate)) / nrOfStates;
	            Arrays.fill(matrix, fPMove);
	            
	            for (int i = 0; i < nrOfStates; i++) {
	                matrix[i * (nrOfStates + 1)] = fPStay;
	            }
	        }
	    }

	    @Override
	    public EigenDecomposition getEigenDecomposition(Node node) {
	        return eigenDecomposition;
	    }

	    @Override
	    public boolean canHandleDataType(DataType dataType) {
	        if (dataType instanceof StandardData || dataType instanceof Binary) {
	            return true;
	        }
	        return false;
	        
	    }

	    protected boolean requiresRecalculation() {
	        if( ! hasFreqs ) {
	            return false;
	        }
	        
	        updateFreqs = true;
	        return true;
	    }
	}
