from Types cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from MSExperiment cimport *
from MzTab cimport *
from StringList cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLModificationsGenerator.h>" namespace "OpenMS":

    cdef cppclass ONuXLModificationsGenerator:

        ONuXLModificationsGenerator() nogil except +
        ONuXLModificationsGenerator(ONuXLModificationsGenerator) nogil except + 

        ONuXLModificationMassesResult initModificationMassesRNA(
            StringList target_nucleotides,
            libcpp_set[char] can_xl,
            StringList mappings,
            StringList modifications,
            String sequence_restriction, 
            bool cysteine_adduct, 
            Int max_length) nogil except +

cdef extern from "<OpenMS/ANALYSIS/XLMS/ONuXLModificationsGenerator.h>" namespace "OpenMS":

    cdef cppclass ONuXLModificationMassesResult:

        ONuXLModificationMassesResult() nogil except +
        ONuXLModificationMassesResult(ONuXLModificationMassesResult) nogil except + 

        # libcpp_map[String, double] mod_masses # empirical formula -> mass
        # libcpp_map[String, libcpp_set[String] ] mod_combinations # empirical formula -> nucleotide formula(s) (formulas if modifications lead to ambiguities)
        # libcpp_map[Size, String] mod_formula_idx


