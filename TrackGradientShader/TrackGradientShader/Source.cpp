#define WIN32_LEAN_AND_MEAN

#include <vector>
#include <fstream>
#include "bmp/EasyBMP.h"
#include "windows.h"
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

const int samples = 1000;
const double thresh = 0.05;

string file1 = "1.txt";
string file2 = "2.txt";
string fileGradient = "gradient.txt";

const int paddingMargin = 10; //padding from track edge to border.

////////////////////////////////////////////////////////////////////////////////

struct Point
{
  int x, y;
  Point(const int& x, const int& y)
    : x(x), y(y) {}
};

////////////////////////////////////////////////////////////////////////////////

BMP canvas;

int max_x = 0;
int max_y = 0;

int percentagetracker = 0;
int totalPointsCount = 0;
double totalPointsProcessed = 0;

////////////////////////////////////////////////////////////////////////////////

bool PopulateVector_Points(const string filename, vector<Point>& v)
{
  ifstream fileio;
  fileio.open(filename);
  
  if (fileio.is_open() == false)
  {
    return false;
  }

  char buff[16] = {0};
  while (!fileio.eof())
  {
    fileio.getline(buff, 16);

    if (strlen(buff) == 0)
    {
      continue;
    }

    char* locator = buff;
    while(*locator!=' ')
    {
      ++locator;
    }
    *locator = '\0';
    ++locator;

    const int x = atoi(buff);
    const int y = atoi(locator);

    v.emplace_back(Point(x, y));

    max_x = (max_x < x) ? x : max_x;
    max_y = (max_y < y) ? y : max_y;
  }

  fileio.close();
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void AdjustPointsForPaddingAndSquareDimensions(vector<Point>& v, const int maxDim)
{
  const int adjust_x = abs(max_x - maxDim)/2 + paddingMargin;
  const int adjust_y = abs(max_y - maxDim)/2 + paddingMargin;

  for (Point& p : v)
  {
    p.x += adjust_x;
    p.y += adjust_y;
  }
}

////////////////////////////////////////////////////////////////////////////////

bool PopulateVector_Gradient(const string filename, vector<int>& v)
{
  ifstream fileio;
  fileio.open(filename);

  if (fileio.is_open() == false)
  {
    return false;
  }

  char buff[16] = { 0 };
  while (!fileio.eof())
  {
    fileio.getline(buff, 16);

    if (strlen(buff) == 0)
    {
      continue;
    }

    v.emplace_back(atoi(buff));
  }

  fileio.close();
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void Paint_Bresenham(const Point& p1, const Point& p2, const RGBApixel& color)
{
  int x0 = p1.x;
  int y0 = p1.y;
  const int x1 = p2.x;
  const int y1 = p2.y;

  const auto dx = abs(x1 - x0);
  const auto dy = abs(y1 - y0);
  const auto sx = (x0 < x1) ? 1 : -1;
  const auto sy = (y0 < y1) ? 1 : -1;
  auto err = dx - dy;

  while (true) 
  {
    canvas.SetPixel(x0, y0, color);

    if ((x0 == x1) && (y0 == y1)) break;
    const auto e2 = 2 * err;
    if (e2 >-dy) { err -= dy; x0 += sx; }
    if (e2 < dx) { err += dx; y0 += sy; }
  }
}

////////////////////////////////////////////////////////////////////////////////

double GetDistance(const Point& p1, const Point& p2)
{
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

////////////////////////////////////////////////////////////////////////////////

void ShadeTrackV1(vector<Point>& v1, vector<Point>& v2, vector<int> gradient)
{
  const int v1Size = v1.size();
  const int gradientLength = gradient.size();

  for (double i = 0; i < v1Size; ++i)
  {
    const int greyval = gradient[(i / v1Size) * gradientLength];
    const RGBApixel colour = { greyval, greyval, greyval, 0 };
    //const RGBApixel colour = { rand() % 255, rand() % 255, rand() % 255, 0 };
    Paint_Bresenham(v1[i], v2[i], colour);
  }
}

////////////////////////////////////////////////////////////////////////////////

void ShadeTrackV2(vector<Point>& v1, vector<Point>& v2, vector<int> gradient)
{
  const int v1Size = v1.size();
  const int v2Size = v2.size();
  const int gradientLength = gradient.size();

  for (double i = 0; i < v1Size; ++i)
  {
    const int greyval = gradient[(i / v1Size) * gradientLength];
    const RGBApixel colour = { greyval, greyval, greyval, 0 };
    //const RGBApixel colour = { rand() % 255, rand() % 255, rand() % 255, 0 };

    const Point p1 = v1[i];
    int p2Index = 0;
    double minDistance = (GetDistance(v1[i], v2[p2Index]));
    for (double j = 0; j < v2Size; ++j)
    {
      if (abs(i - j)>samples)
      {
        continue;
      }

      const Point p2 = v2[j];
      const double dist = GetDistance(p1, p2);
      if (dist < minDistance)
      {
        minDistance = dist < minDistance ? dist : minDistance;
        p2Index = j;
      }
    }

    Paint_Bresenham(p1, v2[p2Index], colour);
  }
}

////////////////////////////////////////////////////////////////////////////////

void ShadeTrackV3(vector<Point>& v1, vector<Point>& v2, vector<int> gradient)
{
  const int v1Size = v1.size();
  const int v2Size = v2.size();
  const int gradientLength = gradient.size();

  for (double i = 0; i < v1Size; ++i)
  {
    const int greyval = gradient[(i / v1Size) * gradientLength];
    const RGBApixel colour = { greyval, greyval, greyval, 0 };
    //const RGBApixel colour = { rand() % 255, rand() % 255, rand() % 255, 0 };

    const Point p1 = v1[i];

    double minDistance = (GetDistance(p1, v2[0]));
    for (double j = 0; j < v2Size; ++j)
    {
      if (abs(i - j)>samples)
      {
        continue;
      }

      const double dist = GetDistance(p1, v2[j]);
      if (dist < minDistance)
      {
        minDistance = dist < minDistance ? dist : minDistance;
      }
    }

    const double minDistanceAdjusted = minDistance + (minDistance * thresh);

    vector<Point> drawablePoints;

    for (auto p : v2)
    {     
      const double distToV2 = GetDistance(p1, p);
      if (distToV2 <= minDistanceAdjusted)
      {
        drawablePoints.push_back(p);
      }
    }

    for (auto pe : drawablePoints)
    {
      Paint_Bresenham(p1, pe, colour);
    }

    const int percentage = (++totalPointsProcessed / totalPointsCount) * 100;
    if (percentage > percentagetracker)
    {
      ++percentagetracker;
      cout << percentage << "%\n";
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

inline bool IsWhite(const RGBApixel& p)
{
  return (p.Red + p.Green + p.Blue) == 765;
}

////////////////////////////////////////////////////////////////////////////////

inline bool AnyAreWhite(vector<RGBApixel>& v)
{
  for (auto p : v)
  {
    if (IsWhite(p))
    {
      return true;
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

void FillInGaps()
{
  for (int x = 1; x < canvas.TellWidth()-1; ++x)
  {
    for (int y = 1; y < canvas.TellHeight()-1; ++y)
    {
      const RGBApixel p = canvas.GetPixel(x, y);
      
      if (IsWhite(p))
      {
        vector<RGBApixel> neighbours = 
        { 
          canvas.GetPixel(x, y - 1), 
          canvas.GetPixel(x, y + 1), 
          canvas.GetPixel(x - 1, y),
          canvas.GetPixel(x + 1, y)
        };

        if (AnyAreWhite(neighbours) == true)
        {
          continue;
        }

        int redSum = 0;
        for (auto n : neighbours)
        {
          redSum += n.Red;
        }

        const int redAvg = redSum / neighbours.size();
        const RGBApixel colour = { redAvg, redAvg, redAvg ,0 };
        canvas.SetPixel(x, y, colour);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void SetGroundLevel()
{
  const RGBApixel black = { 0 };

  for (int x = 0; x < canvas.TellWidth(); ++x)
  {
    for (int y = 0; y < canvas.TellHeight(); ++y)
    {
      const RGBApixel p = canvas.GetPixel(x, y);
      
      if (IsWhite(p))
      {
        canvas.SetPixel(x, y, black);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void ShowError(string fileName)
{
  wstring msg = L"Could not find : ";
  msg += wstring(fileName.begin(), fileName.end());

  MessageBox(
    nullptr,
    msg.c_str(),
    L"Error",
    MB_ICONEXCLAMATION
  );
}

///////////////////////////////////////////////////////////////////////////////

int main(const int argc, char* argv[])
{
  //AL.
  //For debugging!
  #ifdef _DEBUG
  MessageBox(nullptr, L"Attach", L"", 0);
  #endif
  for (int i = 0; i < argc; ++i)
  {
    cout << argv[i] << "\n";
  }
  cout << "\nShading Heightmap...\n";
  //

  if (argc > 1)
  {
    const string trackPath = argv[1];
    file1 = trackPath + "//" + file1;
    file2 = trackPath + "//" + file2;
    const string gradientPath = argv[2];
    fileGradient = gradientPath + "//" + fileGradient;
  }

  vector<Point> v1, v2;
  vector<int> gradient;

  if (PopulateVector_Points(file1, v1) == false)
  {
    ShowError(file1);
    return -1;
  }

  if (PopulateVector_Points(file2, v2) == false)
  {
    ShowError(file2);
    return -1;
  }

  if (PopulateVector_Gradient(fileGradient, gradient) == false)
  {
    ShowError(fileGradient);
    return -1;
  }
  
  totalPointsCount = v1.size() + v2.size();

  const int maxDim = max(max_x, max_y);
  canvas.SetSize(maxDim + (paddingMargin * 2), maxDim + (paddingMargin * 2));

  AdjustPointsForPaddingAndSquareDimensions(v1, maxDim);
  AdjustPointsForPaddingAndSquareDimensions(v2, maxDim);

  ShadeTrackV3(v1, v2, gradient);
  ShadeTrackV3(v2, v1, gradient);

  cout << "\nPLEASE WAIT...\n";

  cout << "\nFILLING IN WHITE GAPS...\n";
  FillInGaps();

  cout << "\nSETTING GROUND LEVEL TO BLACK...\n";
  SetGroundLevel();


  const int tickCount = abs(static_cast<int>(GetTickCount()));
  string filename = to_string(tickCount) + ".bmp";
  cout << "\nSaving to file: " << filename << "\n";
  canvas.WriteToFile(filename.c_str());

  return tickCount;
}

////////////////////////////////////////////////////////////////////////////////
