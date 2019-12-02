#include "headers/Source.h"

#include <iostream>
#include "TApplication.h"

void StandaloneApplication(int argc, char** argv) {

  Source source("a", 5.);
  std::cout << "Hello world\n";

}
// This is the standard "main" of C++ starting
// a ROOT application
int main(int argc, char** argv) {

  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplication(app.Argc(), app.Argv());
  app.Run();

  return 0;
}

