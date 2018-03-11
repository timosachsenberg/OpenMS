// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

/**
  @page TOPP_TOPPView TOPPView

  TOPPView is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzML, mzData, mzXML
  and several other file formats. It also supports viewing data from an %OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC-MS map and a MS raw spectrum:

    @image html TOPPView.png

    More information about TOPPView can be found in the @ref TOPP_tutorial.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_TOPPView.cli
*/

//QT
#include <QtGui/QSplashScreen>
#include <QMessageBox>

//OpenMS
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/SYSTEM/StopWatch.h>

using namespace OpenMS;
using namespace std;

//STL
#include <iostream>
#include <map>
#include <vector>

#include <stdlib.h> 


#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif

class TOPPView :
  public TOPPBase
{
public:
  TOPPView() :
    TOPPBase("TOPPView", "A viewer for mass spectrometry data.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<file>", StringList(), "input file", false, false);
    setValidFormats_("in", ListUtils::create<String>("mzML,featureXML,idXML,consensusXML"));
  }

  ExitCodes main_(int argc, const char ** argv)
  {
    try
    {
      QApplicationTOPP a(argc, const_cast<char**>(argv));
      a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

      TOPPViewBase* mw = new TOPPViewBase();
      a.connect(&a, SIGNAL(fileOpen(QString)), mw, SLOT(loadFile(QString)));
      mw->show();

      // Create the splashscreen that is displayed while the application loads
      QSplashScreen* splash_screen = new QSplashScreen(QPixmap(":/TOPPView_Splashscreen.png"));
      splash_screen->show();
      splash_screen->showMessage("Loading parameters");
      QApplication::processEvents();
      StopWatch stop_watch;
      stop_watch.start();

      //load command line files
      StringList in_files = getStringList_("in");
      if (!in_files.empty())
      {
        mw->loadFiles(in_files, splash_screen);
      }

      // We are about to show the application.
      // Proper time to  remove the splashscreen, if at least 1 second has passed...
      while (stop_watch.getClockTime() < 1.0) {}
      stop_watch.stop();
      splash_screen->close();
      delete splash_screen;

#ifdef OPENMS_WINDOWSPLATFORM
      FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
      AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

      int result = a.exec();
      delete(mw);
      if (!result) return EXECUTION_OK;
      return UNKNOWN_ERROR;
    }
    //######################## ERROR HANDLING #################################
    catch (Exception::UnableToCreateFile& e)
    {
      cout << String("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::FileNotFound& e)
    {
      cout << String("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::FileNotReadable& e)
    {
      cout << String("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::FileEmpty& e)
    {
      cout << String("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::ParseError& e)
    {
      cout << String("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::InvalidValue& e)
    {
      cout << String("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }
    catch (Exception::BaseException& e)
    {
      cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
    }

    return UNKNOWN_ERROR;
  }

};

int main(int argc, const char** argv)
{
  // parse command line into param object and reformat it for CTD support
  Param param;
  Map<String, String> valid_options, valid_flags, option_lists;
  
  valid_options["-ini"] = "ini";
  valid_options["-in"] = "in";
  valid_options["--help"] = "help";
  valid_options["--helphelp"] = "helphelp";
  valid_options["-write_ini"] = "write_ini";
  valid_options["-write_ctd"] = "write_ctd";
  valid_options["-write_wsdl"] = "write_wsdl";
  param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

  if (argc == 1)  // no parameters at all
  {
    TOPPView tool;
    return static_cast<int>(tool.main(argc, argv, true));  // true: ignore empty argument list
  }

  // check for some odd mac parameter
  if (argc > 1 
   && !param.exists("in") 
   && !param.exists("write_ctd") 
   && !param.exists("ini") 
   && !param.exists("help") 
   && !param.exists("helphelp") 
   && !param.exists("write_ini") 
   && !param.exists("write_wsdl"))
  {
    // test if unknown options were given
    if (param.exists("unknown"))
    {
      // if TOPPView is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
      // if this is the only unknown option it will be ignored .. maybe this should be solved directly
      // in Param.h
      if (!(param.getValue("unknown").toString().hasSubstring("-psn") 
       && !param.getValue("unknown").toString().hasSubstring(", ")))
      {
        cout << "Unknown option(s) '" 
             << param.getValue("unknown").toString() 
             << "' given. Aborting!" << endl;
        return 1;
      }
    }

    // make plain command line string compatible with TOPP tool cmd line options
    // prepend -in to command line arguments
    int newc = argc + 1;
    char **newv = (char **)calloc((newc + 1), sizeof(*newv)); // newc + 1 for terminating 0
    memmove(newv + 2, argv + 1, sizeof(*newv) * (argc - 1)); // copy all arguments except the filename to newv[2..newc]
    memmove(newv, argv, sizeof(*newv)); // copy filename (first argument)
    const char* in_string = "-in";
    newv[1] = in_string;
    newv[newc] = 0; // terminating 0

    TOPPView tool;
    return static_cast<int>(tool.main(newc, (const char **)newv));
  }

  // if "-in" or "-ini" parameter is given call the tool with current command line
  if (param.exists("in") 
   || param.exists("ini") 
   || param.exists("help") 
   || param.exists("helphelp")
   || param.exists("write_ini")
   || param.exists("write_ctd"))
  {
    TOPPView tool;
    return static_cast<int>(tool.main(argc, (const char **)argv));  
  }

}

