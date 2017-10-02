//Samuel Wein temporary nucleic acid sequence class


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>
#include <string>


namespace OpenMS
{

NASequence::NASequence()
{ //Default to RNA since that was the original behaviour
    s_ = "";
    type_ = Ribonucleotide::RNA;
}

NASequence::NASequence(const String& rhs)
{ //Default to RNA since that was the original behaviour
    s_ = rhs;
    type_ = Ribonucleotide::RNA;
}

NASequence::NASequence(const String& rhs, const Ribonucleotide::NucleicAcidType& type)
{
    s_ = rhs;
    type_ = type;
}

NASequence& NASequence::operator=(const NASequence& rhs)
{
    s_ = rhs.s_;
    type_ = rhs.type_;
    return *this;
}

bool NASequence::operator==(const NASequence& rhs) const
{
    if (s_ == rhs.getSequence() && type_ == rhs.getType())
        return true;
    else
        return false;
}

NASequence::~NASequence()
{
} //destructor


void NASequence::setSequence(const String& s)
{
    s_ = s;
}

void NASequence::setType(const Ribonucleotide::NucleicAcidType& type)
{
    type_ = type; // @TODO: what if the sequence doesn't fit the type ("T" vs. "U")?
}

Ribonucleotide::NucleicAcidType NASequence::getType() const
{
    return type_;
}

String NASequence::getSequence() const
{
    return s_;
}

bool NASequence::empty() const
{
    if (s_.empty())
        return true;
    else
        return false;
}

size_t NASequence::size() const
{
    size_t size = s_.size();
    // don't count terminal phosphates in the sequence length:
    size_t terminal_p = 0;
    if ((size > 0) && (s_[0] == 'p')) terminal_p++;
    if ((size > 1) && (s_[size - 1] == 'p')) terminal_p++;
    return size - terminal_p;
}

NASequence NASequence::getPrefix(Size index) const
{
    if (index>=s_.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, s_.size());
    return NASequence(s_.substr(0,index),type_);
}

NASequence NASequence::getSuffix(Size index) const
{
    if (index >= s_.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, s_.size());
    return NASequence(s_.substr(s_.size() - index), type_);
}

EmpiricalFormula NASequence::getFormula(Ribonucleotide::RiboNucleotideType type, Int charge) const{

    static const EmpiricalFormula H_weight = EmpiricalFormula("H");
    static const EmpiricalFormula OH_weight = EmpiricalFormula("OH");
    static const EmpiricalFormula NH_weight = EmpiricalFormula("NH");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
    static const EmpiricalFormula fivePrime_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula threePrime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("O"); // OH or just O?
    static const EmpiricalFormula d_ion_to_full = EmpiricalFormula("H2O"); //H2O falls off here
    static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("O");
    static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula w_ion_to_full = EmpiricalFormula("");
    EmpiricalFormula abasicform;
    if (type_ == Ribonucleotide::DNA)
    {
        abasicform = EmpiricalFormula("C5H7O5P");
    }
    else
    {
        abasicform = EmpiricalFormula("C5H7O6P");
    }
    //generate monophosphate mass list TODO make this static
    std::map<char, EmpiricalFormula> base_to_formula;
    //remove H2O since it gets added in internal_to_full
    base_to_formula['A'] = EmpiricalFormula("C5H5N5"); //("C10H12N5O6P");
    base_to_formula['C'] = EmpiricalFormula("C4H5N3O"); //("C9H12N3O7P");
    base_to_formula['B'] = EmpiricalFormula("C5H7N3O"); //
    base_to_formula['G'] = EmpiricalFormula("C5H5N5O"); //("C10H12N5O7P");
    base_to_formula['#'] = EmpiricalFormula("C6H7N5O"); // 2'-O-methyl G
    base_to_formula['T'] = EmpiricalFormula("C5H6N2O2"); //
    base_to_formula['U'] = EmpiricalFormula("C4H4N2O2"); //("C9H11N2O8P");
    base_to_formula['J'] = EmpiricalFormula("C5H6N2O2"); //2'-O-methyl U
    base_to_formula['p'] = EmpiricalFormula("HPO3"); //Placeholder for terminal phosphate
    //C5H7O6P= PO4
    EmpiricalFormula mono_formula;

    // double mono_weight(Constants::PROTON_MASS_U * charge*-1); //the original assumed positive mode

    if (s_.size() > 0)
    {
        if (s_.size() == 0) //FIXME
        {
            mono_formula += base_to_formula[s_[0]] + abasicform;
            return mono_formula + (H_weight * charge) + internal_to_full; //FIXME add switch for other terminals (phosphates etc.)
        }
        else
        {
            for (size_t i = 0; i < s_.size(); i++)
            {
                if (s_[i] == 'p')
                    mono_formula += base_to_formula[s_[i]]; //special case to handle terminal phosphates
                else
                    mono_formula += base_to_formula[s_[i]] + abasicform;
            }
            //            for (ConstIterator it = s_->begin(); it != s_->end(); ++it)
            //            {
            //                // standard residue including named modifications
            //                mono_weight += it->getMonoWeight(Residue::Internal);
            //            }

            // add the missing formula part
            switch (type)
            {
            case Ribonucleotide::Full:
                return mono_formula + internal_to_full - fivePrime_to_full + (H_weight * charge);

                //            case Residue::Internal:
                //                return EmpiricalFormula("");//mono_formula-(H_weight*charge) /* THIS IS NOT CORRECT AND SHOULDNT BE USED FIXME*/;

            case Ribonucleotide::FivePrime:
                return mono_formula + internal_to_full - fivePrime_to_full + (H_weight * charge);

                //            case Residue::CTerminal:
                //                return EmpiricalFormula("");//mono_formula + internal_to_full - threePrime_to_full-(H_weight*charge); //NEED TO CHECK WHAT IS CORRECT FIXME

            case Ribonucleotide::BIon:
                return mono_formula + internal_to_full - b_ion_to_full - H_weight + (H_weight * charge);

            case Ribonucleotide::AIon:
                return mono_formula + internal_to_full - a_ion_to_full - H_weight * 2 + (H_weight * charge);

            case Ribonucleotide::CIon:
                return mono_formula + internal_to_full - c_ion_to_full + (H_weight * charge);

            case Ribonucleotide::DIon:
                return mono_formula + internal_to_full - d_ion_to_full + (H_weight * charge);

            case Ribonucleotide::XIon:
                return mono_formula + internal_to_full - x_ion_to_full + (H_weight * charge);

            case Ribonucleotide::WIon:
                return mono_formula + internal_to_full - w_ion_to_full + (H_weight * charge);

            case Ribonucleotide::YIon:
                return mono_formula + internal_to_full - y_ion_to_full + (H_weight * charge);

            case Ribonucleotide::ZIon:
                return mono_formula + internal_to_full - z_ion_to_full + (H_weight * charge);

            case Ribonucleotide::AminusB:
                return mono_formula + internal_to_full - a_ion_to_full - H_weight * 2 + (H_weight * charge) - base_to_formula[s_[s_.size()-1]];

            default:
                LOG_ERROR << "NASequence::getMonoWeight: unknown RibonucleotideType" << std::endl;
            }
        }
    }

    return mono_formula;
}

double NASequence::getMonoWeight(Ribonucleotide::RiboNucleotideType type, Int charge) const
{
    //    static const double H_weight = EmpiricalFormula("H").getMonoWeight();
    //    static const double OH_weight = EmpiricalFormula("OH").getMonoWeight();
    //    static const double NH_weight = EmpiricalFormula("NH").getMonoWeight();
    //    static const double internal_to_full = EmpiricalFormula("H2O").getMonoWeight();
    //    static const double fivePrime_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double threePrime_to_full = EmpiricalFormula("H").getMonoWeight();
    //    static const double b_ion_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double a_ion_to_full = EmpiricalFormula("PO4").getMonoWeight();
    //    static const double c_ion_to_full = EmpiricalFormula("O").getMonoWeight(); // OH or just O?
    //    static const double d_ion_to_full = EmpiricalFormula("H2O").getMonoWeight(); //H2O falls off here
    //    static const double x_ion_to_full = EmpiricalFormula("O").getMonoWeight();
    //    static const double y_ion_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double z_ion_to_full = EmpiricalFormula("PO4").getMonoWeight();
    //    static const double w_ion_to_full = EmpiricalFormula("").getMonoWeight();
    //    //generate monophosphate mass list TODO make this static
    //    std::map<char, EmpiricalFormula> monophosphate_to_formula;
    //    //remove H2O since it gets added in internal_to_full
    //    monophosphate_to_formula['A']=EmpiricalFormula("C10H12N5O6P");
    //    monophosphate_to_formula['C']=EmpiricalFormula("C9H12N3O7P");
    //    monophosphate_to_formula['G']=EmpiricalFormula("C10H12N5O7P");
    //    monophosphate_to_formula['U']=EmpiricalFormula("C9H11N2O8P");

    double mono_weight(getFormula(type, charge).getMonoWeight()); //the original assumed positive mode
    return mono_weight;//(double)charge;//+getFormula(type,charge).getMonoWeight();

    //    if (s_.size() > 0)
    //    {
    //        if (s_.size() == 1)
    //        {
    //            return mono_weight + monophosphate_to_formula[s_[0]].getMonoWeight();
    //        }
    //        else
    //        {
    //            for (size_t i=0;i<s_.size();i++){
    //                mono_weight+=monophosphate_to_formula[s_[i]].getMonoWeight();
    //            }
    ////            for (ConstIterator it = s_->begin(); it != s_->end(); ++it)
    ////            {
    ////                // standard residue including named modifications
    ////                mono_weight += it->getMonoWeight(Residue::Internal);
    ////            }

    //            // add the missing formula part
    //            switch (type)
    //            {
    //            case Residue::Full:
    //                return mono_weight + internal_to_full;

    //            case Residue::Internal:
    //                return mono_weight /* + add_protons*/;

    //            case Residue::NTerminal:
    //                return mono_weight + internal_to_full - fivePrime_to_full;

    //            case Residue::CTerminal:
    //                return mono_weight + internal_to_full - threePrime_to_full;

    //            case Residue::BIon:
    //                return mono_weight + internal_to_full - b_ion_to_full - H_weight;

    //            case Residue::AIon:
    //                return mono_weight + internal_to_full - a_ion_to_full - H_weight * 2;

    //            case Residue::CIon:
    //                return mono_weight + internal_to_full - c_ion_to_full;

    //            case Residue::DIon:
    //                return mono_weight + internal_to_full - d_ion_to_full;

    //            case Residue::XIon:
    //                return mono_weight + internal_to_full - x_ion_to_full;

    //            case Residue::WIon:
    //                return mono_weight + internal_to_full - w_ion_to_full;

    //            case Residue::YIon:
    //                return mono_weight + internal_to_full - y_ion_to_full;

    //            case Residue::ZIon:
    //                return mono_weight + internal_to_full - z_ion_to_full;

    //            default:
    //                LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
    //            }
    //        }
    //    }


    //    return mono_weight;
}

OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide){
    os << nucleotide.getSequence();
    return os;
}
}