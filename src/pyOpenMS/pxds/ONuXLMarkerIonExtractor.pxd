from Types cimport *
from MSSpectrum cimport *
from String cimport *

# typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;

ctypedef libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] MarkerIonsType

cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLMarkerIonExtractor.h>" namespace "OpenMS":
    
    cdef cppclass ONuXLMarkerIonExtractor "OpenMS::ONuXLMarkerIonExtractor":
        ONuXLMarkerIonExtractor() nogil except + 
        ONuXLMarkerIonExtractor(ONuXLMarkerIonExtractor) nogil except + #wrap-ignore
        # MarkerIonsType extractMarkerIons(MSSpectrum & s, double marker_tolerance) nogil except +
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] extractMarkerIons(MSSpectrum & s,
                                                                                            double marker_tolerance) nogil except +

