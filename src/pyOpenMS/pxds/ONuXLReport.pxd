from Types cimport *

from ONuXLMarkerIonExtractor cimport *
# from ListUtilsIO cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass ONuXLReport "OpenMS::ONuXLReport":
        ONuXLReport(ONuXLReport) nogil except + #wrap-ignore
        libcpp_vector[ ONuXLReportRow ] annotate(MSExperiment & spectra,
                                                 libcpp_vector[ PeptideIdentification ] & peptide_ids,
                                                 double marker_ions_tolerance) nogil except +


cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass ONuXLReportRowHeader "OpenMS::ONuXLReportRowHeader":
        ONuXLReportRowHeader(ONuXLReportRowHeader) nogil except + #wrap-ignore
        String getString(const String & separator) nogil except +

cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass ONuXLReportRow "OpenMS::ONuXLReportRow":
        ONuXLReportRow(ONuXLReportRow) nogil except + #wrap-ignore
        bool no_id
        double rt
        double original_mz
        String accessions
        String RNA
        String peptide
        double best_localization_score
        String localization_scores
        String best_localization
        Int charge
        double score
        double peptide_weight
        double RNA_weight
        double xl_weight
        double abs_prec_error
        double rel_prec_error
        # ONuXLMarkerIonExtractor::MarkerIonsType marker_ions
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] marker_ions
        double m_H
        double m_2H
        double m_3H
        double m_4H
        String getString(const String & separator) nogil except +

