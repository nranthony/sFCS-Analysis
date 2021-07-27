#include<algorithm>
#include <vector>
#include "mex.h"


/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    /* Check for proper number of arguments */
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("MATLAB:sfcs_carpet:nargin", "sfcs_carpet requires seven input arguments: sfcsCarpet(uint32(deltat), uint32(Nc), uint32(C), uint16(microt), uint64(macrot), uint64(marktime), uint8(marktype))");
    }

    uint32_t* pdt = (uint32_t*) mxGetData(prhs[0]);
    uint32_t* pNc = (uint32_t*) mxGetData(prhs[1]);
    uint32_t* pC  = (uint32_t*) mxGetData(prhs[2]);

    uint32_t dt = *pdt;
    uint32_t Nc = *pNc;
    uint32_t C = *pC;

    uint16_t* pmicrot = (uint16_t*) mxGetData(prhs[3]);
    uint64_t* pmacrot = (uint64_t*) mxGetData(prhs[4]);
    uint64_t* pmrkrstime = (uint64_t*) mxGetData(prhs[5]);
    uint8_t* pmrkrstype = (uint8_t*) mxGetData(prhs[6]);

    std::vector<std::vector<uint64_t>> macrotVects(Nc);
    std::vector<std::vector<uint16_t>> microtVects(Nc);

    // inc marker pointers to first '4' - assumes sequence 1, 4, 2, 4, ... , 2, 4
    pmrkrstime++;
    pmrkrstype++;

    uint8_t mrktype = NULL;
    uint64_t mrktime = NULL;
    uint64_t pixmrktime = NULL;
    uint16_t mict = NULL;
    uint64_t mact = NULL;
    int64_t delta_mact = NULL;

    plhs[0] = mxCreateNumericMatrix((mwSize)C, (mwSize)Nc, mxUINT8_CLASS, mxREAL);
    uint8_t* pcounts = (uint8_t*) mxGetData(plhs[0]);
    //uint8_t* pcounts = new uint8_t[C * Nc];
    uint8_t cnt = 0;
    bool in_time_window = true;

    // inc macro pointer until inside the first time window
    while (*pmacrot < *pmrkrstime) {
        pmacrot++;
        pmicrot++;
    }

    // TODO - if *pmacrot < mrktime things never progress...

    uint32_t inWindowCount = 0;

    for (size_t c = 0; c < C; c++)
    {   
        // TODO - consider sfcs vs linefcs with these markers
        mrktime = *(pmrkrstime + (2 * c)); // set temp marker at start of each cycle in C

        for (size_t n = 0; n < Nc; n++)
        {
            in_time_window = true;
            cnt = 0; // reset counter for inserting into counts
            pixmrktime = mrktime + (n * dt); // offset for each 'pixel' in Nc cycle points
            delta_mact = (int64_t)*pmacrot - (int64_t)pixmrktime; 
            if (delta_mact < 0) {
                // sometimes a fixed delta t for each pixel and timing jitter on markers can lead to tiny gaps
                // assign photon to previous window, inc pointers and continue
                // debugging:
                //mexPrintf("*pmacrot (%u) - pixmrktime (%u) is negative (%i)\n", *pmacrot, pixmrktime, delta_mact);
                //return;
                cnt++;
                inWindowCount++;
                size_t w = NULL;
                if (n == 0) { w = Nc - 1; }
                else { w = n - 1; }
                macrotVects[w].push_back(*pmacrot);
                microtVects[w].push_back(*pmicrot);
                pmacrot++;
                pmicrot++;
            }
            while (in_time_window)
            {
                // think about reversing the process: get macrot, check window, inc until it fits, inc macro
                if ((*pmacrot >= pixmrktime) & (*pmacrot < (pixmrktime + dt))) {
                    // macrot between marker time and deltat
                    cnt++;
                    inWindowCount++;
                    macrotVects[n].push_back(*pmacrot);
                    microtVects[n].push_back(*pmicrot);
                    pmacrot++;
                    pmicrot++;
                }
                else {
                    // macrot ahead of window - move to next; record count
                    in_time_window = false;
                    // pcounts ptr to uint8_t[C * Nc]
                    *(pcounts + ((c * Nc) + n)) = cnt;
                }
            }
        }
    }

    mexPrintf("%u photon counts in sliding windows.\nLast pixel window %i - %i\n", inWindowCount, pixmrktime, pixmrktime + dt);

    mxArray* cell_array_ptr = mxCreateCellMatrix((mwSize)2, (mwSize)Nc);
    size_t cell_counter = 0;

    for (size_t n = 0; n < Nc; n++)
    {
        mxArray* macroArr = mxCreateNumericMatrix(macrotVects[n].size(), 1, mxUINT64_CLASS, mxREAL);
        uint64_t* pmacArr = (uint64_t*) mxGetData(macroArr);
        uint64_t* pmacVect = macrotVects[n].data();
 
        mxArray* microArr = mxCreateNumericMatrix(microtVects[n].size(), 1, mxUINT16_CLASS, mxREAL);
        uint16_t* pmicArr = (uint16_t*)mxGetData(microArr);
        uint16_t* pmicVect = microtVects[n].data();
        
        //mexPrintf("\n")

        for (size_t i = 0; i < macrotVects[n].size(); i++)
        {   
            pmacArr[i] = pmacVect[i];
            pmicArr[i] = pmicVect[i];
        }
        mxSetCell(cell_array_ptr, cell_counter, macroArr);
        cell_counter++;
        mxSetCell(cell_array_ptr, cell_counter, microArr);
        cell_counter++;
        
    }

    plhs[1] = cell_array_ptr;

}
