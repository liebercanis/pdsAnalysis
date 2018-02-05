#include "display.hxx"
display* disp;
void eventDisplay()
{
  disp = new display();
  //disp->AddTag("07-31-2059_0");
  //disp->AddTag("07-31-1614_0");
  disp->AddTag("07-31-1646_0");
  //disp->AddTag("07-31-1714_0");
  //disp->AddTag("07-31-2020_0");
  //disp->AddTag("07-31-1518_0");
  //disp->AddTag("07-31-1700_0");
  disp->SetTagList();
  printf(" type disp->DrawEvent(number)\n");
 //   disp->DrawEvent(8);
 //   disp->print(".tiff");
 //   disp->DrawEvent(18);
 //   disp->print(".tiff");
 return;
  
  for(int i=1066; i<1069; ++i) {
    disp->DrawEvent(i);
    disp->print(".tiff");
  }
}
