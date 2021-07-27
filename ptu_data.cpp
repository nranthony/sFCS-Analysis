#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

//#define header_debug

#pragma pack(8) //structure alignment to 8 byte boundaries

// some important Tag Idents (TTagHead.Ident) that we will need to read the most common content of a PTU file
// check the output of this program and consult the tag dictionary if you need more

#define TTTRTagTTTRRecType      "TTResultFormat_TTTRRecType"    
#define TTTRTagNumRecords       "TTResult_NumberOfRecords"      // Number of TTTR Records in the File
#define TTTRTagAcqTime          "MeasDesc_AcquisitionTime"      // Acquistion Time in ms
#define TTTRTagPixX             "ImgHdr_PixX"                   // Number of pixels in 'x' direction - see labels
#define TTTRTagPixY             "ImgHdr_PixY"                   // Number of pixels in 'y' direction - see labels
#define TTTRTagTimePerPix       "ImgHdr_TimePerPixel"           // Time per pixel in ms
#define TTTRTagScanAxLabel0     "UsrAIScanAx_Label(0)"          // Measurement axis (0) label
#define TTTRTagScanAxLabel1     "UsrAIScanAx_Label(1)"          // Measurement axis (1) label
#define TTTRTagPixRes           "ImgHdr_PixResol"               // Physicsal pixel resolution
#define TTTRTagSyncRate         "TTResult_SyncRate"             // Laser pulse sync rate frequency in Hz
#define TTTRTagBiDirect         "ImgHdr_BiDirect"               // Boolean for bidirectional scanning
#define TTTRTagRes              "MeasDesc_Resolution"           // Resolution for the Dtime (T3 Only)
#define TTTRTagGlobRes          "MeasDesc_GlobalResolution"     // Global Resolution of TimeTag(T2) /NSync (T3)
#define FileTagEnd              "Header_End"                    // Always appended as last tag (BLOCKEND)

// TagTypes  (TTagHead.Typ)
#define tyEmpty8      0xFFFF0008
#define tyBool8       0x00000008
#define tyInt8        0x10000008
#define tyBitSet64    0x11000008
#define tyColor8      0x12000008
#define tyFloat8      0x20000008
#define tyTDateTime   0x21000008
#define tyFloat8Array 0x2001FFFF
#define tyAnsiString  0x4001FFFF
#define tyWideString  0x4002FFFF
#define tyBinaryBlob  0xFFFFFFFF

// TDateTime (in file) to time_t (standard C) conversion

const int EpochDiff = 25569; // days between 30/12/1899 and 01/01/1970
const int SecsInDay = 86400; // number of seconds in a day

time_t TDateTime_TimeT(double Convertee)
{
    time_t Result((long)(((Convertee)-EpochDiff) * SecsInDay));
    return Result;
}

using namespace matlab::data;
using matlab::mex::ArgumentList;
using std::cout;
using std::endl;
using std::hex;


class MexFunction : public matlab::mex::Function {
public:

    FILE* fpin;
    char Magic[8];
    char Version[8];
    char Buffer[40];
    char* AnsiBuffer;
    wchar_t* WideBuffer;
    int Result;
    bool PrintHeader;

    double GlobRes = 0.0;
    double Resolution = 0.0;
    unsigned int dlen = 0;
    unsigned int cnt_0 = 0, cnt_1 = 0;

    // header parameters
    int RecordLength = 4; // Length of an Record, by default 4 bytes, But LIN Camera has 8 bytes
    long long RecordType = 0;
    long long NumRecords = -1;
    unsigned long AcqTime = 0;
    unsigned short PixX = 0;
    unsigned short PixY = 0;
    double TimePerPix = 0.0;
    std::string ScanAxis0 = "";
    std::string ScanAxis1 = "";
    double PixRes = 0.0;
    unsigned long SyncRate = 0;
    bool BiDirect = false;

    ~MexFunction() {
        if (fpin != NULL)
        {
           fclose(fpin);
        }        
    }


    void operator()(ArgumentList outputs, ArgumentList inputs) {

        const double inScalar = inputs[1][0];
        bool prnthdr = (bool)inScalar;
        PrintHeader = prnthdr;

        ArrayFactory af;
        
        std::ostringstream stream;

        // TODO - check large arrays to test if this is stack limited
        //  using 'new' defines vector on the heap
        std::vector<uint64_t> macrotVect;
        //std::vector<uint64_t>* macrotVect = new std::vector<uint64_t>();
        std::vector<uint16_t> microtVect;
        //std::vector<uint16_t>* microtVect = new std::vector<uint16_t>();

        std::vector<uint64_t> markerTimeVect;
        std::vector<uint8_t> markerTypeVect;

        // Validate arguments
        // TODO - check filepath exists here
        //checkArguments(outputs, inputs);

        // Read header and populate parameters
        readHeader(inputs, stream);

        // Pull TTTR records in vectors
        ProcessTTTRData(macrotVect, microtVect, markerTimeVect, markerTypeVect);

        // true ns per pulse - leave as units of sync/laser pulses for now 
        //double tns;
        //tns = 1000000000.0 / (double)SyncRate;

        // these objects are std::unique_ptr<T[],matlab::data::buffer_deleter_t>
        auto Vmac_b_ptr = af.createBuffer<uint64_t>(macrotVect.size());
        // this pointer points to the first adress that the std::unique_ptr points to
        uint64_t* Vmac_ptr = Vmac_b_ptr.get();
        // multiply each element by tns to give truetime in ns - leave for now - units actually number of laser pulses
        //std::transform(macrotVect.begin(), macrotVect.end(), macrotVect.begin(), [&](uint64_t v) { return std::round( (double)v * tns); });
        // fill the buffer with the V values
        std::for_each(macrotVect.begin(), macrotVect.end(), [&](const uint64_t& e) { *(Vmac_ptr++) = e; });
        // finally create Matlab array from buffer
        TypedArray<uint64_t> Amac = af.createArrayFromBuffer<uint64_t>({ macrotVect.size(), 1 }, std::move(Vmac_b_ptr));
        outputs[0] = Amac;

        // these objects are std::unique_ptr<T[],matlab::data::buffer_deleter_t>
        auto Vmic_b_ptr = af.createBuffer<uint16_t>(microtVect.size());
        // this pointer points to the first adress that the std::unique_ptr points to
        uint16_t* Vmic_ptr = Vmic_b_ptr.get();
        // fill the buffer with the V values
        std::for_each(microtVect.begin(), microtVect.end(), [&](const uint16_t& e) { *(Vmic_ptr++) = e; });
        // finally create Matlab array from buffer
        TypedArray<uint16_t> Amic = af.createArrayFromBuffer<uint16_t>({ microtVect.size(), 1 }, std::move(Vmic_b_ptr));
        outputs[1] = Amic;

        // histogram microtime data
        auto dec_b_ptr = af.createBuffer<uint32_t>(1024);
        uint32_t* dec_ptr = dec_b_ptr.get();
        for (size_t i = 0; i < 1024; i++) *(dec_ptr + i) = 0;
        uint32_t n;
        std::for_each(microtVect.begin(), microtVect.end(), [&](const uint16_t& e) {
            n = *(dec_ptr + e);
            n++;
            *(dec_ptr + e) = n;
            });
        TypedArray<uint32_t> Adec = af.createArrayFromBuffer<uint32_t>({ 1024, 1 }, std::move(dec_b_ptr));
        outputs[2] = Adec;

        // TODO: 
        // create and fill a buffer for both marker time and type - insert into outputs[3] and outputs[4] respectively
        // 
        // markerTimeVect -> matlab buffer
        auto Vmktime_b_ptr = af.createBuffer<uint64_t>(markerTimeVect.size());
        // this pointer points to the first adress that the std::unique_ptr points to
        uint64_t* Vmktime_ptr = Vmktime_b_ptr.get();
        // fill the buffer with the V values
        std::for_each(markerTimeVect.begin(), markerTimeVect.end(), [&](const uint64_t& e) { *(Vmktime_ptr++) = e; });
        // finally create Matlab array from buffer
        TypedArray<uint64_t> Amktime = af.createArrayFromBuffer<uint64_t>({ markerTimeVect.size(), 1 }, std::move(Vmktime_b_ptr));
        outputs[3] = Amktime;

        // markerTypeVect -> matlab buffer
        auto Vmktype_b_ptr = af.createBuffer<uint8_t>(markerTypeVect.size());
        // this pointer points to the first adress that the std::unique_ptr points to
        uint8_t* Vmktype_ptr = Vmktype_b_ptr.get();
        // fill the buffer with the V values
        std::for_each(markerTypeVect.begin(), markerTypeVect.end(), [&](const uint8_t& e) { *(Vmktype_ptr++) = e; });
        // finally create Matlab array from buffer
        TypedArray<uint8_t> Amktype = af.createArrayFromBuffer<uint8_t>({ markerTypeVect.size(), 1 }, std::move(Vmktype_b_ptr));
        outputs[4] = Amktype;


        // add variables to cell array
        if (outputs.size() == 6) {
            //  optional parameter cell array requested
            CellArray metaCellArr = af.createCellArray({ 2,10 }
                //TTTRTagNumRecords   long long NumRecords = -1;    // Number of TTTR Records in the File
                , af.createScalar<int64_t>(NumRecords), af.createScalar("NumRecords")
                //TTTRTagAcqTime      unsigned long AcqTime = 0;    // Acquistion Time in ms
                , af.createScalar<uint32_t>(AcqTime), af.createScalar("AcqTime(ms)")
                //TTTRTagPixX         unsigned short PixX = 0;      // Number of pixels in 'x' direction
                , af.createScalar<uint16_t>(PixX), af.createScalar("PixX")
                //TTTRTagPixY         unsigned short PixY = 0;      // Number of pixels in 'y' direction
                , af.createScalar<uint16_t>(PixY), af.createScalar("PixY")
                //TTTRTagTimePerPix   double TimePerPix = 0.0;      // Time per pixel in ms  
                , af.createScalar<double>(TimePerPix), af.createScalar("TimePerPix(ms)")
                //TTTRTagPixRes       double PixRes = 0.0;          // Physicsal pixel resolution
                , af.createScalar<double>(PixRes), af.createScalar("PixRes(um)")
                //TTTRTagSyncRate     unsigned long SyncRate = 0;   // Laser pulse sync rate frequency in Hz
                , af.createScalar<uint32_t>(SyncRate), af.createScalar("SyncRate(Hz)")
                //TTTRTagBiDirect     bool BiDirect = false;        // Boolean for bidirectional scanning
                , af.createScalar<bool>(BiDirect), af.createScalar("BiDirect")
                //TTTRTagRes          double Resolution = 0.0;      // Resolution for the Dtime (T3 Only)
                , af.createScalar<double>(Resolution), af.createScalar("TCSPC_ChnRes(s)")
                //TTTRTagGlobRes      double GlobRes = 0.0;         // Global Resolution of TimeTag(T2)/NSync (T3)
                , af.createScalar<double>(GlobRes), af.createScalar("Pule2Pulse(s)")     
            );                                        
            outputs[5] = metaCellArr;
        }

    }

    void readHeader(ArgumentList inputs, std::ostringstream& stream) {
        const CharArray fpath = inputs[0];
        std::string fpath_ascii;
        fpath_ascii = fpath.toAscii();

        char outStr[128];
        int numOutStr;

        // A Tag entry
        struct TgHd {
            char Ident[32];     // Identifier of the tag
            int Idx;            // Index for multiple tags or -1
            unsigned int Typ;  // Type of tag ty..... see const section
            long long TagValue; // Value of tag.
        } TagHead;

        stream << fpath_ascii << std::endl;
        displayOnMATLAB(stream);

        if ((fpin = fopen(fpath_ascii.c_str(), "rb")) == NULL)
        {
            numOutStr = sprintf(outStr, "ERROR! Input file cannot be opened, aborting.");
            stream << outStr << std::endl;
            displayOnMATLAB(stream);
            return;
        }

        // Check ptu file
        Result = fread(&Magic, 1, sizeof(Magic), fpin);
        if (Result != sizeof(Magic))
        {
            numOutStr = sprintf(outStr, "Error reading header, aborted.");
            stream << outStr << std::endl;
            displayOnMATLAB(stream);
            return;
        }
        Result = fread(&Version, 1, sizeof(Version), fpin);
        if (Result != sizeof(Version))
        {
            numOutStr = sprintf(outStr, "Error reading header, aborted.");
            stream << outStr << std::endl;
            displayOnMATLAB(stream);
            return;
        }
        if (strncmp(Magic, "PQTTTR", 6))
        {
            numOutStr = sprintf(outStr, "Wrong Magic, this is not a PTU file.");
            stream << outStr << std::endl;
            displayOnMATLAB(stream);
            return;
        }
        numOutStr = sprintf(outStr, "PQTTTR File, Tag Version: %s ", Version);
        stream << outStr << std::endl;
        displayOnMATLAB(stream);

        // read tagged header
        do
        {
            // This loop is very generic. It reads all header items and displays the identifier and the
            // associated value, quite independent of what they mean in detail.
            // Only some selected items are explicitly retrieved and kept in memory because they are
            // needed to subsequently interpret the TTTR record data.

            // Need to check for these and populate the variables from header
            // TODO - return parameters in struct or cell array
            // using ~ to indicate completed tags; not sure how to handle the array types - leave for now
            /*
            ~TTTRTagTTTRRecType      long long RecordType = 0;     //
            ~TTTRTagNumRecords       long long NumRecords = -1;    // Number of TTTR Records in the File
            ~TTTRTagAcqTime          unsigned long AcqTime = 0;    // Acquistion Time in ms
            ~TTTRTagPixX             unsigned short PixX = 0;      // Number of pixels in 'x' direction - see labels
            ~TTTRTagPixY             unsigned short PixY = 0;      // Number of pixels in 'y' direction - see labels
            ~TTTRTagTimePerPix       double TimePerPix = 0.0;      // Time per pixel in ms
            TTTRTagScanAxLabel0     string ScanAxis0 = "";        // Measurement axis (0) label
            TTTRTagScanAxLabel1     string ScanAxis1 = "";        // Measurement axis (1) label
            ~TTTRTagPixRes           double PixRes = 0.0;          // Physicsal pixel resolution
            ~TTTRTagSyncRate         unsigned long SyncRate = 0;   // Laser pulse sync rate frequency in Hz
            ~TTTRTagBiDirect         bool BiDirect = false;        // Boolean for bidirectional scanning
            ~TTTRTagRes              double Resolution = 0.0;         // Resolution for the Dtime (T3 Only)
            ~TTTRTagGlobRes          double GlobRes = 0.0;      // Global Resolution of TimeTag(T2) /NSync (T3)
            */
            Result = fread(&TagHead, 1, sizeof(TagHead), fpin);
            if (Result != sizeof(TagHead))
            {
                numOutStr = sprintf(outStr, "Incomplete File.");
                stream << outStr << std::endl;
                displayOnMATLAB(stream);
                return;
            }

            strcpy(Buffer, TagHead.Ident);
            if (TagHead.Idx > -1)
            {
                numOutStr = sprintf(Buffer, "%s(%d)", TagHead.Ident, TagHead.Idx);
            }
            numOutStr = sprintf(outStr, "%-40s", Buffer);
            stream << outStr;
            switch (TagHead.Typ)
            {
            case tyEmpty8:
                numOutStr = sprintf(outStr, "<empty Tag>");
#ifdef header_debug
                stream << outStr << "[tyEmpty8]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                break;
            case tyBool8:
                numOutStr = sprintf(outStr, "%s", bool(TagHead.TagValue) ? "True" : "False");              
#ifdef header_debug
                stream << outStr << "[tyBool8]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                if (strcmp(TagHead.Ident, TTTRTagBiDirect) == 0) // Boolean for bidirectional scanning
                    BiDirect = *(bool*)&(TagHead.TagValue);
                break;
            case tyInt8:
                numOutStr = sprintf(outStr, "%lld", TagHead.TagValue);            
#ifdef header_debug
                stream << outStr << "[tyInt8]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                // get some Values we need to analyse records
                if (strcmp(TagHead.Ident, TTTRTagNumRecords) == 0) // Number of records
                    NumRecords = TagHead.TagValue;
                if (strcmp(TagHead.Ident, TTTRTagTTTRRecType) == 0) // TTTR RecordType
                    RecordType = TagHead.TagValue;
                if (strcmp(TagHead.Ident, TTTRTagAcqTime) == 0) //  Acq timein ms
                    AcqTime = TagHead.TagValue;
                if (strcmp(TagHead.Ident, TTTRTagPixX) == 0) //  Number of pixels in 'x' direction - see labels
                    PixX = TagHead.TagValue;
                if (strcmp(TagHead.Ident, TTTRTagPixY) == 0) //  Number of pixels in 'y' direction - see labels
                    PixY = TagHead.TagValue;
                if (strcmp(TagHead.Ident, TTTRTagSyncRate) == 0) // Laser pulse sync rate frequency in Hz
                    SyncRate = TagHead.TagValue;
                break;
            case tyBitSet64:
                numOutStr = sprintf(outStr, "0x%16.16X", TagHead.TagValue);                
#ifdef header_debug
                stream << outStr << "[tyBitSet64]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                break;
            case tyColor8:
                numOutStr = sprintf(outStr, "0x%16.16X", TagHead.TagValue);                
#ifdef header_debug
                stream << outStr << "[tyColor8]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                break;
            case tyFloat8:
                numOutStr = sprintf(outStr, "%E", *(double*)&(TagHead.TagValue));                
#ifdef header_debug
                stream << outStr << "[tyFloat8]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                if (strcmp(TagHead.Ident, TTTRTagRes) == 0) // Resolution for TCSPC-Decay
                    Resolution = *(double*)&(TagHead.TagValue);
                if (strcmp(TagHead.Ident, TTTRTagGlobRes) == 0) // Global resolution for timetag
                    GlobRes = *(double*)&(TagHead.TagValue); // in ns
                if (strcmp(TagHead.Ident, TTTRTagTimePerPix) == 0) // Time per pixel in ms
                    TimePerPix = *(double*)&(TagHead.TagValue);
                if (strcmp(TagHead.Ident, TTTRTagPixRes) == 0) // Physicsal pixel resolution
                    PixRes = *(double*)&(TagHead.TagValue);
                break;
            case tyFloat8Array:
                numOutStr = sprintf(outStr, "<Float Array with %d Entries>", TagHead.TagValue / sizeof(double));                
#ifdef header_debug
                stream << outStr << "[tyFloat8Array]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                // only seek the Data, if one needs the data, it can be loaded here
                fseek(fpin, (long)TagHead.TagValue, SEEK_CUR);
                break;
            case tyTDateTime:
                time_t CreateTime;
                CreateTime = TDateTime_TimeT(*((double*)&(TagHead.TagValue)));
                numOutStr = sprintf(outStr, "%s", asctime(gmtime(&CreateTime)), "\0");               
#ifdef header_debug
                stream << outStr << "[tyTDateTime]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                break;
            case tyAnsiString:
                AnsiBuffer = (char*)calloc((size_t)TagHead.TagValue, 1);
                Result = fread(AnsiBuffer, 1, (size_t)TagHead.TagValue, fpin);
                if (Result != TagHead.TagValue)
                {
                    numOutStr = sprintf(outStr, "Incomplete File.");
                    stream << outStr << std::endl;
                    displayOnMATLAB(stream);
                    free(AnsiBuffer);
                    return;
                }
                numOutStr = sprintf(outStr, "%s", AnsiBuffer);               
#ifdef header_debug
                stream << outStr << "[tyAnsiString]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                free(AnsiBuffer);
                break;
            case tyWideString:
                WideBuffer = (wchar_t*)calloc((size_t)TagHead.TagValue, 1);
                Result = fread(WideBuffer, 1, (size_t)TagHead.TagValue, fpin);
                if (Result != TagHead.TagValue)
                {
                    numOutStr = sprintf(outStr, "Incomplete File.");                    
#ifdef header_debug
                    stream << outStr << "[tyWideString]" << std::endl;
#else
                    stream << outStr << std::endl;
#endif
                    displayOnMATLAB(stream);
                    free(WideBuffer);
                    return;
                }
                numOutStr = sprintf(outStr, "%ws", WideBuffer);
                stream << outStr << std::endl;
                free(WideBuffer);
                break;
            case tyBinaryBlob:
                numOutStr = sprintf(outStr, "<Binary Blob contains %d Bytes>", TagHead.TagValue);                
#ifdef header_debug
                stream << outStr << "[tyBinaryBlob]" << std::endl;
#else
                stream << outStr << std::endl;
#endif
                // only seek the Data, if one needs the data, it can be loaded here
                fseek(fpin, (long)TagHead.TagValue, SEEK_CUR);
                break;
            default:
                numOutStr = sprintf(outStr, "Illegal Type identifier found! Broken file?");
                stream << outStr << std::endl;
                displayOnMATLAB(stream);
                return;
            }

            displayOnMATLAB(stream);

        }         while ((strncmp(TagHead.Ident, FileTagEnd, sizeof(FileTagEnd))));
    }

    void ProcessTTTRData(
        std::vector<uint64_t> &macrotVect,
        std::vector<uint16_t> &microtVect,
        std::vector<uint64_t> &markerTimeVect,
        std::vector<uint8_t> &markerTypeVect )
    {
        
        unsigned int TTTRRecord;
        const int T3WRAPAROUND = 1024;
        long long oflcorrection = 0;
        long long RecNum;
        long long truensync, truetime;
        unsigned short m;
        uint8_t c;


        union {
            unsigned int allbits;
            struct {
                unsigned nsync : 10;    // numer of sync period
                unsigned dtime : 15;    // delay from last sync in units of chosen resolution
                unsigned channel : 6;
                unsigned special : 1;
            } bits;
        } T3Rec;
        T3Rec.allbits = TTTRRecord;

        for (RecNum = 0; RecNum < NumRecords; RecNum++)
        {
            Result = fread(&TTTRRecord, 1, sizeof(TTTRRecord), fpin);
            T3Rec.allbits = TTTRRecord;

            if (T3Rec.bits.special == 1)
            {
                if (T3Rec.bits.channel == 0x3F) //overflow
                {
                    oflcorrection += (int64_t)T3WRAPAROUND * T3Rec.bits.nsync;
                }
                else if ((T3Rec.bits.channel >= 1) && (T3Rec.bits.channel <= 15)) //markers
                {
                    truensync = oflcorrection + T3Rec.bits.nsync;
                    markerTimeVect.push_back(truensync);
                    c = T3Rec.bits.channel;
                    markerTypeVect.push_back(c);

                }
            }
            else //regular input channel
            {
                truensync = oflcorrection + T3Rec.bits.nsync;
                macrotVect.push_back(truensync);

                //the nsync time unit depends on sync period which can be obtained from the file header
                //the dtime unit depends on the resolution and can also be obtained from the file header
                //c = T3Rec.bits.channel;  // assuming single channel for the moment
                m = T3Rec.bits.dtime;
                microtVect.push_back(m);
            }
        }
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory af;
        if (inputs[0].getType() != ArrayType::DOUBLE || inputs[0].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error", 0, std::vector<Array>({ af.createScalar("Input must be double array") }));
        }
        if (outputs.size() > 1)
        {
            matlabPtr->feval(u"error", 0, std::vector<Array>({ af.createScalar("Only one output is returned") }));
        }
    }

    void displayOnMATLAB(std::ostringstream& stream) {
        if (PrintHeader) {
            std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
            ArrayFactory af;
            // Pass stream content to MATLAB fprintf function
            matlabPtr->feval("fprintf", 0, std::vector<Array>({ af.createScalar(stream.str()) }));
            // Clear stream buffer
            stream.str("");
        }
    }
    
};