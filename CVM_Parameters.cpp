/******   Scenario   ******/
constexpr static char Scenario = 'B'; //  Scenario A/B/C  

/******   Outputfiles   ******/
constexpr char DATAOUTFILE[] = "results/result_data.txt" ; // 
constexpr int  GNUPLOT_USE = 2;  //  0: no gnuplot  1: show gnuplot 2: save as gif animation
constexpr char OUTGIFFILE[]  = "results/cvm.gif";  // gif 

/******   Main variables; Loaded from initial configuration file   ******/
int CELL_NUMBER; 
int V_NUM; // (2*CELL_NUMBER)
int E_NUM; // (3*CELL_NUMBER)
int N;  // (2*V_NUM)

vector<double> X; // vertices positions 
double *x, *y; // pointer 
double LX, LY; // system size, set by loading initial configuraiton file


/******   CVM parameters   ******/
double BVo;    // BV0 initalized in the main function A0 = L^2/N
double KV  = 1.0;  // volume elasticity, 
double KL  = 1.0e-1; // 1.0e-2; // 0.01; //
double ShpIdx = 3.7; //  shape index  hexagon 3.72 
double BLo;          // L0 is set as ShpIdx * (A0)^{1/2}
vector<double> Vo, Lo; // Vo and Lo for individual cells.

/******   Cell division parameters   ******/
// We suppose 1 time unit in this simulation = 2 min.
// 'm' cells divide and change to two 'd' cells after MITOTIC_TIME x2 min.
constexpr double MITOTIC_TIME = 10.0; // mitotic phase duration, = 20 min. 
int cell_div_counter = 0, div_by_nucleation = 0, div_by_triggered = 0; // # of transition to 'm', that by nucleation and by neigbour triggering


/***  Simulation time and output intertval ***/
constexpr int SIMULATION_TIME = 250;
constexpr int TIME_INTERVAL =  1000;
