
#include "mex.h"
#include <string>

static std::string ss = "disp('";

namespace Cantera {

    void writelog(const std::string& s) {
        char ch = s[0];
        int n = 0;
        while (ch != '\0') {
            if (ch =='\n') {
                ss += "');";
                mexEvalString(ss.c_str());
                ss = "disp('";
            }
            else 
                ss += ch;
            if (ch == '\'') ss += ch;
            n++;
            ch = s[n];
        }
    }
    
}
