// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"

// STV analysis includes
#include "DAgostiniUnfolder.hh"
#include "FiducialVolume.hh"
#include "MCC9SystematicsCalculator.hh"
#include "NormShapeCovMatrix.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"
#include "WienerSVDUnfolder.hh"
#include "Constants.hh"

using namespace Constants;

//----------------------------------------//

void multiply_1d_hist_by_matrix( TMatrixD* mat, TH1* hist ) {

  // Copy the histogram contents into a column vector

  int num_bins = mat->GetNcols();
  TMatrixD hist_mat( num_bins, 1 );

  for ( int r = 0; r < num_bins; ++r ) {

    hist_mat( r, 0 ) = hist->GetBinContent( r + 1 );

  }

  // Multiply the column vector by the input matrix
  // TODO: add error handling here related to matrix dimensions
  TMatrixD hist_mat_transformed( *mat, TMatrixD::EMatrixCreatorsOp2::kMult,hist_mat );

  // Update the input histogram contents with the new values
  for ( int r = 0; r < num_bins; ++r ) {

    double val = hist_mat_transformed( r, 0 );
    hist->SetBinContent( r + 1, val );

  }

}

//----------------------------------------//

struct SampleInfo {

  SampleInfo() {}

  SampleInfo( const std::string& resp, const std::string& width,
    const std::string& sb ) : respmat_file_( resp ), widths_file_( width ),
    sb_file_( sb ) {}

  // File containing the output of respmat.C with the same true binning
  // as was used in NUISANCE to make theoretical predictions
  std::string respmat_file_;

  // NUISANCE data file in which the bin widths (along each relevant axis) are
  // stored. These files are used by the preliminary NUISANCE implementation of
  // this cross-section measurement.
  std::string widths_file_;

  // SliceBinning configuration file (needs to be consistent with the bin
  // definitions used in the other two files)
  std::string sb_file_;
};

//----------------------------------------//

// Keys are NUISANCE sample names, values are SampleInfo objects that give
// the corresponding input file paths needed for this script
std::map< std::string, SampleInfo > sample_info_map = {

  // 2D proton measurement
  { "MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC", { "/uboone/data/users/apapadop/PeLEE_PostProcessing/new_stv_analysis.root", "/uboone/app/users/apapadop/master_stv_analysis/mybins_mcc9_2D_proton.txt",
    "/uboone/app/users/apapadop/master_stv_analysis/myconfig_mcc9_2D_proton.txt" } },

};

//----------------------------------------//

void dump_slice_errors( const std::string& hist_col_prefix,
  const Slice& slice, const std::map< std::string,
  std::unique_ptr<SliceHistogram> >& slice_hist_cov_matrix_map,
  std::map< std::string, std::vector<double> >& pgf_plots_hist_table )
{

  for ( const auto& pair : slice_hist_cov_matrix_map ) {

    std::string err_name = pair.first;
    std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";
    pgf_plots_hist_table[ err_col_name ] = std::vector<double>();

  }

  for ( const auto& bin_pair : slice.bin_map_ ) {

    int global_bin_idx = bin_pair.first;

    for ( const auto& err_pair : slice_hist_cov_matrix_map ) {

      std::string err_name = err_pair.first;
      std::string err_col_name = hist_col_prefix + '_' + err_name + "_error";

      const auto* hist = err_pair.second->hist_.get();
      double err = hist->GetBinError( global_bin_idx );

      pgf_plots_hist_table.at( err_col_name ).push_back( err );
    }

  } // slice bins

}

//----------------------------------------//

void xsec_unfolding() {

  //----------------------------------------//

  gStyle->SetEndErrorSize(6);

  //----------------------------------------//

  // Root file with all the xsec results

  TFile* xsec_file = new TFile(PATH_XSEC+"xsec_analyzer_results.root","recreate");

  //----------------------------------------//

  //// Initialize the FilePropertiesManager and tell it to treat the NuWro
  //// MC ntuples as if they were data
  //auto& fpm = FilePropertiesManager::Instance();
  //fpm.load_file_properties( std::string(WORK_AREA) + "nuwro_file_properties.txt" );

  const auto& sample_info = sample_info_map.at( SAMPLE_NAME );
  const std::string respmat_file_name( PATH_POSTPROCESS + "new_stv_analysis.root" );

  // Do the systematics calculations in preparation for unfolding
  auto* syst_ptr = new MCC9SystematicsCalculator( respmat_file_name, std::string(WORK_AREA) + "systcalc.conf" );
  auto& syst = *syst_ptr;

  // Get the tuned GENIE CV prediction in each true bin (including the
  // background true bins)
  TH1D* genie_cv_truth = syst.cv_universe().hist_true_.get();
  int num_true_bins = genie_cv_truth->GetNbinsX();

  // While we're at it, clone the histogram and zero it out. We'll fill this
  // one with our unfolded result for easy comparison
  TH1D* unfolded_events = dynamic_cast< TH1D* >(
    genie_cv_truth->Clone("unfolded_events") );
  unfolded_events->Reset();

  // If present, then get the fake data event counts in each true bin
  // (including the background true bins). We hope to approximately reproduce
  // these event counts in the signal true bins via unfolding the fake data.
  const auto& fake_data_univ = syst.fake_data_universe();
  TH1D* fake_data_truth_hist = nullptr;

  bool using_fake_data = false;

  if ( fake_data_univ ) {

    using_fake_data = true;
    fake_data_truth_hist = fake_data_univ->hist_true_.get();

  }

  int num_ordinary_reco_bins = 0;
  int num_sideband_reco_bins = 0;

  for ( int b = 0; b < syst.reco_bins_.size(); ++b ) {

    const auto& rbin = syst.reco_bins_.at( b );
    if ( rbin.type_ == kSidebandRecoBin ) ++num_sideband_reco_bins;
    else ++num_ordinary_reco_bins;

  }

  int num_true_signal_bins = 0;
  for ( int t = 0; t < syst.true_bins_.size(); ++t ) {

    const auto& tbin = syst.true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;

  }

  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();

  constexpr int NUM_DAGOSTINI_ITERATIONS = 2;
  constexpr bool USE_ADD_SMEAR = true;

  std::unique_ptr< Unfolder > unfolder (
    // //new DAgostiniUnfolder( NUM_DAGOSTINI_ITERATIONS )
    //new DAgostiniUnfolder( DAgostiniUnfolder::ConvergenceCriterion
    //  ::FigureOfMerit, 0.025 )
    new WienerSVDUnfolder( true,
      WienerSVDUnfolder::RegularizationMatrixType::kSecondDeriv )
  );

  UnfoldedMeasurement result = unfolder->unfold( syst );

  // For real data only, add some new covariance matrices in which only the
  // signal response or the background is varied. We could calculate these
  // for the fake data, but it seems unnecessary at this point.
  if ( !using_fake_data ) {

    syst.set_syst_mode( MCC9SystematicsCalculator
      ::SystMode::VaryOnlyBackground );
    auto* bkgd_matrix_map_ptr = syst.get_covariances().release();
    auto& bkgd_matrix_map = *bkgd_matrix_map_ptr;

    syst.set_syst_mode( MCC9SystematicsCalculator
      ::SystMode::VaryOnlySignalResponse );
    auto* sigresp_matrix_map_ptr = syst.get_covariances().release();
    auto& sigresp_matrix_map = *sigresp_matrix_map_ptr;

    for ( const auto& m_pair : bkgd_matrix_map ) {

      auto& my_temp_cov_mat = matrix_map[ "bkgd_only_" + m_pair.first ];
      my_temp_cov_mat += m_pair.second;

    }

    for ( const auto& m_pair : sigresp_matrix_map ) {

      auto& my_temp_cov_mat = matrix_map[ "sigresp_only_" + m_pair.first ];
      my_temp_cov_mat += m_pair.second;

    }

  }

  // Propagate all defined covariance matrices through the unfolding procedure
  const TMatrixD& err_prop = *result.err_prop_matrix_;
  TMatrixD err_prop_tr( TMatrixD::kTransposed, err_prop );

  std::map< std::string, std::unique_ptr<TMatrixD> > unfolded_cov_matrix_map;

  for ( const auto& matrix_pair : matrix_map ) {

    const std::string& matrix_key = matrix_pair.first;
    auto temp_cov_mat = matrix_pair.second.get_matrix();

    TMatrixD temp_mat( *temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult,
      err_prop_tr );

    unfolded_cov_matrix_map[ matrix_key ] = std::make_unique< TMatrixD >(
      err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

  }

  // Decompose the block-diagonal pieces of the total covariance matrix
  // into normalization, shape, and mixed components (for later plotting
  // purposes)
  NormShapeCovMatrix bd_ns_covmat = make_block_diagonal_norm_shape_covmat(
    *result.unfolded_signal_, *result.cov_matrix_, syst.true_bins_ );

  // Add the blockwise decomposed matrices into the map
  unfolded_cov_matrix_map[ "total_blockwise_norm" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.norm_ );

  unfolded_cov_matrix_map[ "total_blockwise_shape" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.shape_ );

  unfolded_cov_matrix_map[ "total_blockwise_mixed" ]
    = std::make_unique< TMatrixD >( bd_ns_covmat.mixed_ );

  for (auto const &pair: unfolded_cov_matrix_map) {
    std::cout << pair.first << endl;
  }

  // Set the event counts in each bin of the histogram that displays the
  // unfolded result. Note that we don't care about the background true bins
  // (which are assumed to follow all of the signal true bins) since we've
  // subtracted out an estimate of the background before unfolding.

  for ( int t = 0; t < num_true_bins; ++t ) {

    double evts = 0.;
    double error = 0.;

    if ( t < num_true_signal_bins ) {

      evts = result.unfolded_signal_->operator()( t, 0 );
      error = std::sqrt( std::max(0., result.cov_matrix_->operator()( t, t )) );

    }

    // We need to use one-based indices while working with TH1D bins
    unfolded_events->SetBinContent( t + 1, evts );
    unfolded_events->SetBinError( t + 1, error );

  }

  unfolded_events->SetStats( false );
  unfolded_events->SetLineColor( kBlack );
  unfolded_events->SetLineWidth( LINE_THICKNESS );
  unfolded_events->GetXaxis()->SetRangeUser( 0, num_true_signal_bins );

  // Multiply the truth-level GENIE prediction by the additional smearing
  // matrix
  TMatrixD* A_C = result.add_smear_matrix_.get();
  multiply_1d_hist_by_matrix( A_C, genie_cv_truth );

  if ( using_fake_data ) {

    // Multiply the fake data truth by the additional smearing matrix
    multiply_1d_hist_by_matrix( A_C, fake_data_truth_hist );

  }

  // Plot slices of the unfolded result
  auto* sb_ptr = new SliceBinning( std::string(WORK_AREA) + "tutorial_slice_config.txt" );
  auto& sb = *sb_ptr;

  // Get the factors needed to convert to cross-section units
  double total_pot = syst.total_bnb_data_pot_;
  double integ_flux = integrated_numu_flux_in_FV( total_pot );
  double num_Ar = num_Ar_targets_in_FV();

  // Retrieve the true-space expected event counts from NUISANCE output files
  // for each available generator model
  double conv_factor = num_Ar * integ_flux;

  // Add the fake data truth using a column vector of event counts
  TMatrixD fake_data_truth( num_true_signal_bins, 1 );

  if ( using_fake_data ) {

    for ( int b = 0; b < num_true_signal_bins; ++b ) {

      double true_evts = fake_data_truth_hist->GetBinContent( b + 1 );
      fake_data_truth( b, 0 ) = true_evts;

    }

  }

  // Add the GENIE CV model using a column vector of event counts
  TMatrixD genie_cv_truth_vec( num_true_signal_bins, 1 );
  for ( int b = 0; b < num_true_signal_bins; ++b ) {

    double true_evts = genie_cv_truth->GetBinContent( b + 1 );
    genie_cv_truth_vec( b, 0 ) = true_evts;

  }

  if ( USE_ADD_SMEAR ) {

    // Get access to the additional smearing matrix
    const TMatrixD& A_C = *result.add_smear_matrix_;

    // Start with the fake data truth if present
    if ( using_fake_data ) {

      TMatrixD ac_truth( A_C, TMatrixD::kMult, fake_data_truth );
      fake_data_truth = ac_truth;

    }

    // Also transform the GENIE CV model
    TMatrixD genie_cv_temp( A_C, TMatrixD::kMult, genie_cv_truth_vec );
    genie_cv_truth_vec = genie_cv_temp;

  }

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // Make a histogram showing the unfolded true event counts in the current
    // slice
    SliceHistogram* slice_unf = SliceHistogram::make_slice_histogram(
      *result.unfolded_signal_, slice, result.cov_matrix_.get() );

    // Temporary copies of the unfolded true event count slices with
    // different covariance matrices
    std::map< std::string, std::unique_ptr<SliceHistogram> > sh_cov_map;
    for ( const auto& uc_pair : unfolded_cov_matrix_map ) {

      const auto& uc_name = uc_pair.first;
      const auto& uc_matrix = uc_pair.second;

      auto& uc_ptr = sh_cov_map[ uc_name ];
      uc_ptr.reset(
        SliceHistogram::make_slice_histogram( *result.unfolded_signal_, slice,
        uc_matrix.get() )
      );

    }

    // Also use the GENIE CV model to do the same
    SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(
      genie_cv_truth_vec, slice, nullptr );

    // If present, also use the truth information from the fake data to do the
    // same
    SliceHistogram* slice_truth = nullptr;
    if ( using_fake_data ) {

      slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth,
        slice, nullptr );

    }

    // Keys are legend labels, values are SliceHistogram objects containing
    // true-space predictions from the corresponding generator models
    auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
    auto& slice_gen_map = *slice_gen_map_ptr;

    slice_gen_map[ "Data" ] = slice_unf;
    if ( using_fake_data ) {
      slice_gen_map[ "truth" ] = slice_truth;
    }
    slice_gen_map[ "G18" ] = slice_cv;

    int var_count = 0;
    std::string diff_xsec_denom;
    std::string diff_xsec_units_denom;
    std::string diff_xsec_denom_latex;
    std::string diff_xsec_units_denom_latex;
    double other_var_width = 1.;

    for ( const auto& ov_spec : slice.other_vars_ ) {

      double high = ov_spec.high_bin_edge_;
      double low = ov_spec.low_bin_edge_;
      const auto& var_spec = sb.slice_vars_.at( ov_spec.var_index_ );

      if ( high != low && std::abs(high - low) < BIG_DOUBLE ) {

        ++var_count;
        other_var_width *= ( high - low );
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;
        const std::string& temp_units = var_spec.units_;

        if ( !temp_units.empty() ) {

          diff_xsec_units_denom += " / " + temp_units;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;

        }

      }

    }

    for ( size_t av_idx : slice.active_var_indices_ ) {

      const auto& var_spec = sb.slice_vars_.at( av_idx );
      const std::string& temp_name = var_spec.name_;

      if ( temp_name != "true bin number" ) {

        var_count += slice.active_var_indices_.size();
        diff_xsec_denom += 'd' + var_spec.name_;
        diff_xsec_denom_latex += " d" + var_spec.latex_name_;

        if ( !var_spec.units_.empty() ) {

          diff_xsec_units_denom += " / " + var_spec.units_;
          diff_xsec_units_denom_latex += " / " + var_spec.latex_units_;

        }

      }

    }

    // NOTE: This currently assumes that each slice is a 1D histogram
    // TODO: revisit as needed
    int num_slice_bins = slice_unf->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );

    for ( int b = 0; b < num_slice_bins; ++b ) {

      double width = slice_unf->hist_->GetBinWidth( b + 1 );
      width *= other_var_width;
      trans_mat( b, b ) = 1e38 / ( width * integ_flux * num_Ar );

    }

    std::string slice_y_title;
    std::string slice_y_latex_title;

    if ( var_count > 0 ) {

      slice_y_title += "d";
      slice_y_latex_title += "{$d";

      if ( var_count > 1 ) {

        slice_y_title += "^{" + std::to_string( var_count ) + "}";
        slice_y_latex_title += "^{" + std::to_string( var_count ) + "}";

      }

      slice_y_title += "#sigma/" + diff_xsec_denom;
      slice_y_latex_title += "\\sigma / " + diff_xsec_denom_latex;

    }
    else {

      slice_y_title += "#sigma";
      slice_y_latex_title += "\\sigma";

    }

    slice_y_title += " (10^{-38} cm^{2}" + diff_xsec_units_denom + " / Ar)";
    slice_y_latex_title += "\\text{ }(10^{-38}\\text{ cm}^{2}"
      + diff_xsec_units_denom_latex + " / \\mathrm{Ar})$}";

    // Convert all slice histograms from true event counts to differential
    // cross-section units
    for ( auto& pair : slice_gen_map ) {

      auto* slice_h = pair.second;
      slice_h->transform( trans_mat );
      slice_h->hist_->GetYaxis()->SetTitle( slice_y_title.c_str() );

    }

    // Also transform all of the unfolded data slice histograms which have
    // specific covariance matrices
    for ( auto& sh_cov_pair : sh_cov_map ) {

      auto& slice_h = sh_cov_pair.second;
      slice_h->transform( trans_mat );

    }

    // Keys are generator legend labels, values are the results of a chi^2
    // test compared to the unfolded data (or, in the case of the unfolded
    // data, to the fake data truth)
    std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
    //std::cout << '\n';

    for ( const auto& pair : slice_gen_map ) {

      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      // Decide what other slice histogram should be compared to this one,
      // then calculate chi^2
      SliceHistogram* other = nullptr;
      // We don't need to compare the unfolded data to itself, so just skip to
      // the next SliceHistogram and leave a dummy Chi2Result object in the map
      if ( name == "Data" ) {

        chi2_map[ name ] = SliceHistogram::Chi2Result();
        continue;

      }
      // Compare all other distributions to the unfolded data
      else {

        other = slice_gen_map.at( "Data" );

      }

      // Store the chi^2 results in the map
      const auto& chi2_result = chi2_map[ name ] = slice_h->get_chi2( *other );

      //std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
      //  << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
      //if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
      //std::cout << ", p-value = " << chi2_result.p_value_ << '\n';

    }

    TCanvas* c1 = new TCanvas;
    c1->SetTitle("");
    c1->SetBottomMargin( 0.14 );
    c1->SetLeftMargin( 0.14 );

    slice_unf->hist_->SetLineColor( kBlack );
    slice_unf->hist_->SetLineWidth( LINE_THICKNESS );
    slice_unf->hist_->SetMarkerStyle( kFullCircle );
    slice_unf->hist_->SetMarkerSize( 1. );
    slice_unf->hist_->SetStats( false );
    slice_unf->hist_->SetTitle("");

    slice_unf->hist_->GetXaxis()->CenterTitle();
    slice_unf->hist_->GetXaxis()->SetTitleSize(TEXT_SIZE);
    slice_unf->hist_->GetXaxis()->SetLabelSize(TEXT_SIZE);
    slice_unf->hist_->GetXaxis()->SetTitleFont(TEXT_FONT);
    slice_unf->hist_->GetXaxis()->SetLabelFont(TEXT_FONT);
    slice_unf->hist_->GetXaxis()->SetNdivisions(NDIVISIONS);

    slice_unf->hist_->GetYaxis()->CenterTitle();
    slice_unf->hist_->GetYaxis()->SetTitleSize(TEXT_SIZE);
    slice_unf->hist_->GetYaxis()->SetLabelSize(TEXT_SIZE);
    slice_unf->hist_->GetYaxis()->SetTitleFont(TEXT_FONT);
    slice_unf->hist_->GetYaxis()->SetLabelFont(TEXT_FONT);
    slice_unf->hist_->GetYaxis()->SetNdivisions(NDIVISIONS);

    slice_unf->hist_->GetYaxis()->SetRangeUser(0.,1.15*slice_unf->hist_->GetMaximum());
    slice_unf->hist_->Draw( "e1x0" );

    slice_cv->hist_->SetStats( false );
    slice_cv->hist_->SetLineColor( kAzure - 4 );
    slice_cv->hist_->SetLineWidth( LINE_THICKNESS );
    slice_cv->hist_->SetLineStyle( kSolid );
    slice_cv->hist_->Draw( "hist same" );

    if ( using_fake_data ) {

      slice_truth->hist_->SetStats( false );
      slice_truth->hist_->SetLineColor( kOrange );
      slice_truth->hist_->SetLineWidth( LINE_THICKNESS );
      slice_truth->hist_->Draw( "hist same" );

    }

    slice_unf->hist_->Draw( "e1x0 same" );

    //----------------------------------------//

    // Chi2 & p-value calculation

    TLegend* lg = new TLegend( 0.12, 0.91, 0.9, 0.99 );

    for ( const auto& pair : slice_gen_map ) {

      const auto& name = pair.first;
      const auto* slice_h = pair.second;

      std::string label = name;
      std::ostringstream oss;
      const auto& chi2_result = chi2_map.at( name );
      oss << std::setprecision( 3 ) << chi2_result.chi2_ << "/"
	  << chi2_result.num_bins_ << ", p = " << chi2_result.p_value_ ;

      if ( name != "Data" ) {

        //label += ": #chi^{2} = " + oss.str();
        label += " (#chi^{2}/bins = " + oss.str() + ")";

      }

      if ( name != "Data" ) { lg->AddEntry( slice_h->hist_.get(), label.c_str(), "l" ); }
      else { lg->AddEntry( slice_h->hist_.get(), label.c_str(), "ep" );  }

    }

    lg->SetMargin(0.1);
    lg->SetNColumns(2);
    lg->SetFillColor(0);
    lg->SetBorderSize(0);
    lg->SetTextSize(TEXT_SIZE);
    lg->SetTextFont(TEXT_FONT);
    lg->Draw( "same" );

    //c1->SaveAs(PATH_XSEC + TString( slice_unf->hist_->GetName() ) + ".pdf" );

  } // slices

  return;

}

//----------------------------------------//

int main() {

  xsec_unfolding();
  return 0;

}

//----------------------------------------//
