package simpleMyelomaABM;
import HAL.GridsAndAgents.*;
import HAL.Gui.*;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import static simpleMyelomaABM.simpleBoneGrid.*;
import static HAL.Util.*;
import java.util.concurrent.*;
import java.util.concurrent.ConcurrentLinkedQueue;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     GRID CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public class simpleBoneGrid extends AgentGrid2D<simpleBoneCell> implements SerializableModel {

    ///////////////
    //GRID FIELDS//
    ///////////////


    //Variables to switch on/off treatment/other things
    public static boolean BORTEZOMIB = false;

    public static boolean MYELOMA = true;
    public static boolean TCELL = true;

    public static boolean ADAPTIVE = false;

    public static boolean BIWEEKLY = false;

    public static boolean CIRCULATING = false;
    public static boolean EMDR = false;

    public static boolean TREATMENT_ON = false; //this is to control treatment on/off timer in MAIN

    //CLUSTER
    public static boolean PARAM_SWEEP = false; //use when importing parameters to loop through
    public static boolean runPar = true;

    public static boolean runPareShare = false;
    public static boolean HEADLESS = false; //use true with cluster
    public static boolean LOCAL = true; // use false with cluster
    public static double numSteps = 2*365.0*24.0*60.0; // years the model will run
    public static int numSims = 10; //Number of Simulations
    public final static int BONE = RGB256(255,255,250), MSC = RGB256(135,206,250),
            pOB = RGB256(100,149,237), aOB = BLUE, pOC = RGB256(230,100,130),
            aOC = RED, LINING = RGB256(64,106,151), MM = RGB256(0,128,0),
            activeTcell = RGB256(17, 150, 150),
            EXHT_CELL=RGB256(200, 50, 250),
            supressorTcell =RGB256(255, 165, 0),
            bloodVessel=RGB256(138, 3, 3),
            naiveTcell=RGB256(100, 100, 100);

    //SETUP
    static double MinToHour = 60.0;
    public final static double SPACESTEP = 10.0;//um
    public static double TIMESTEP_AGENT = 6.0/MinToHour; //0.1;//hr; //6.0 min; 6.0/60.0 hour
    public final static double N_TIMESTEP_PDE = 60.0*(MinToHour*TIMESTEP_AGENT);//360.0; //Number of diffusion timesteps; 1 dts = 1 sec; 360 dts = 1 ts = 6 min
    public static int timeStep = 0;

    //DiffCoef MUST <0.25 for FTCS scheme!

    //CHEMOTAXIS
    double Tcell_DiffCoef = 0.01*3.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);
    double Tcell_TaxisCoeff = 5.0e10*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);

    //CELL PARAMETERS
    public int TURNOVER_TIME = (int) (2102400.0/(MinToHour*TIMESTEP_AGENT)); //2102400 min = 350400 ts = 4 years
    public double MM_DEATH = 1.0 / 11000 * (MinToHour);
    public double MM_DEATH_BDF = 1.0 / 72000 * (MinToHour);
    public double pmutate = 0.0;
    public double antigenLoss = Math.pow(10, -2);
    public double T_CELL_DIV_RATE = 1.0 / 720 * (MinToHour); // T_CELL DIVISION RATE
    public double EXHT_CELL_DIV_RATE = (1.0 / 1440)*2 * (MinToHour); // T_CELL DIVISION RATE

    public double TREG_CELL_DIV_RATE = 1.0 / 1440 * (MinToHour); // T_CELL DIVISION RATE

    //Cytokines
    double TGFB_productionRate = 2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//2.04e-9//dts; RANKL production rate
    static double TGFB_basalRate = 2.04e-11*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //2.04e-11//dts; basal TGFB production rate (not by OC)
    static double TGFB_decayRate = -0.35*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //dts; BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    double TGFB_DiffCoef = 780.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE); //dts
    double Extra_TGFBtime = 4320.0/(MinToHour*TIMESTEP_AGENT); //4320 min = 720 ts = 3 days
    static double maxTGFB = 8.7e-10; //1.12e-9; //1.78e-9;//1.21e-9; //6.4e-10

    double Ts = TGFB_basalRate/(Math.abs(TGFB_decayRate)*maxTGFB); //Basal TGFB
    double TGFBthresh = (1.05)* TGFB_basalRate/(Math.abs(TGFB_decayRate)* maxTGFB);//was 0.01; now 5% increase from basal

    //CHEMOTAXIS

    double CXCL9_productionRate = (2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE)/8;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE; //changed from 2.61e-10
    double CXCL9_decayRate = -.2*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE;
    double CXCL9_DiffCoef = 2700.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);//*(TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);
    double maxCXCL9 = (2.07e-9)/8;

    public static int daysPassed = 0;
    //MODEL TESTS
    public double MarrowArea;
    double convert_to_days = (MinToHour*TIMESTEP_AGENT)/(60.0*24.0); //1 ts = 6 min = 1/240 day
    int count_BA = 0;
    int init_BA = 0;
    int Nts = (int) ((numSteps)/(MinToHour*TIMESTEP_AGENT));

    double MaxTdiff; //max relative change in TGFB
    double tmax;//=0;
    public int[] moveHood = VonNeumannHood(true); //Have option of no movement

    public Rand rn;
    public PDEGrid2D TGFB;
    public PDEGrid2D CXCL9;
    public PDEGrid2D IFNG;

    public int[] tmoveHood = MooreHood(true);
    public ArrayList<Integer> InitBoneList = new ArrayList<>();
    public ArrayList<simpleBoneCell> AllBoneList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<simpleBoneCell> LiningList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs

    boolean [][] extraTGFB = new boolean[xDim][yDim];
    double [][] TGFBtimer = new double[xDim][yDim];

    double[] TGFBvals = new double[xDim*yDim];//new double[xDim*yDim];

    List<Object> cellSIMindex = new ArrayList<>();

    //public int eOpt;
    FileIO out;

    FileIO clones;

    FileIO bones;

    FileIO locations;

    FileIO TGFBout;

    FileIO paramsOut;

    FileIO InitialBone;
//    FileIO params;

    //This is important for serializable model
    @Override
    public void SetupConstructors(){
        this._PassAgentConstructor(simpleBoneCell.class);
    }

    ////////////////////
    //GRID CONSTRUCTOR//
    ////////////////////

    public simpleBoneGrid(int xDim, int yDim, Rand rn, String Bone_FileName) {
        super(xDim, yDim, simpleBoneCell.class,true,true);
        this.rn = rn;

        //Create 2D PDE Grids
        TGFB = new PDEGrid2D(xDim, yDim,true,true);
        CXCL9 = new PDEGrid2D(xDim, yDim,true,true);
        IFNG = new PDEGrid2D(xDim, yDim,true,true);
        InitialBone=new FileIO(Bone_FileName, "r");

    }

    /////////////////////////
    //GRID METHODS///////////
    /////////////////////////
    //1. InitBone          //
    //2. RemodelingEvent   //
    //3. InitRANKL         //
    //4. ModelStep         //
    //5. RecordOut         //
    /////////////////////////

    // taking a snapshot of the cell
    public static class CellSnapshot {
        public final int type;
        public final boolean bcmaLoss;
        public int simID;

        public CellSnapshot(int type,
                            boolean bcmaLoss, int simID) {
            this.type = type;
            this.bcmaLoss = bcmaLoss;
            this.simID = simID;
        }
    }

    // --- Add snapshot export helper in BoneGrid_2022May17 (instance method) ---
    public CellSnapshot exportSnapshot(simpleBoneCell c) {
        // Capture the fields you care about. Extend as needed.
        return new CellSnapshot(
                c.type,
                c.bcmaLoss,
                c.simulationID
        );
    }

    // --- Add snapshot import helper in BoneGrid_2022May17 (instance method) ---
    public void importSnapshotToGrid(CellSnapshot s) {
        // Find a random empty location to place the incoming cell (or use other placement logic)
        int tries = 0;
        while (tries < 50) { // try up to 50 random spots
            int x = rn.Int(xDim);
            int y = rn.Int(yDim);
            if (PopAt(x, y) == 0) {
                simpleBoneCell newCell = NewAgentSQ(x, y);
                newCell.type = s.type;
                newCell.bcmaLoss = s.bcmaLoss;
                newCell.simulationID = s.simID;
                // copy other fields as needed
                return;
            }
            tries++;
        }
        // If not placed after random tries, place it at first empty by scanning
        for (int xi = 0; xi < xDim; xi++) {
            for (int yi = 0; yi < yDim; yi++) {
                if (PopAt(xi, yi) == 0) {
                    simpleBoneCell newCell = NewAgentSQ(xi, yi);
                    newCell.type = s.type;
                    newCell.bcmaLoss = s.bcmaLoss;
                    newCell.simulationID = s.simID;
                    return;
                }
            }
        }
    }

    // global, accessible from all threads
    public static ConcurrentLinkedQueue<CellSnapshot> myelomaTransferQueue = new ConcurrentLinkedQueue<>();

    public void IterateTGFB(Grid2Ddouble TGFB_xDiffArray, Grid2Ddouble TGFB_yDiffArray) {
        //Diffusion of TGFB
        TGFB.Diffusion(TGFB_xDiffArray,TGFB_yDiffArray);


        //Natural Decay of TGFB
        TGFB.MulAll(TGFB_decayRate);

        //Production of TGFB
        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                TGFB.Add(x, y, TGFB_basalRate / maxTGFB);
                if (GetAgent(x, y) != null && GetAgent(x, y).TGFB_on == true) {
                    TGFB.Add(x, y, TGFB_productionRate / maxTGFB);
                } else if (extraTGFB[x][y] == true && TGFBtimer[x][y] < Extra_TGFBtime && GetAgent(x, y) == null) {
                    TGFB.Add(x, y, 0.1 * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                }
            }
        }


        //Solved using Euler's method--need small timestep so that concentration doesn't go negative
        if (TGFB.GetMin() < 0) {
            System.out.println("WARNING: Negative TGFB concentration");
        }

        MaxTdiff = TGFB.MaxDeltaScaled(1e-18);
        TGFB.Update(); //This step is necessary to update diffusion each time-step

        //Record TGFB output

    }


    public void GridTick(){
        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if (extraTGFB[x][y] == true ) {
                    TGFBtimer[x][y]++;
                }
                if(TGFBtimer[x][y]>=Extra_TGFBtime && GetAgent(x,y)==null){
                    extraTGFB[x][y]=false;
                    TGFBtimer[x][y]=0;
                }
            }
        }
    }


    public static void runSimulation(final int simID, final int prow, final String baseFolder, final ArrayList<String> param_list, int[] condition) {

        // local dims (same as before)
        final int xDim = 160;
        final int yDim = 150;

        // timestamped subfolder per sim (same naming you used)
        String fn = baseFolder; // base folder created in main
        String subfolder = fn + "/Sim" + simID + "_row" + prow + "/";
        File dir = new File(subfolder);
        dir.mkdirs();

        // -------------------------
        // UI / visualization setup
        // -------------------------
        UIWindow win = HEADLESS ? null : new UIWindow("Normal Bone Remodeling");
        // Full domain UIGrids (same as your deleted code)
        UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
        UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 2);
        UIGrid TGFB_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
        UILabel days = new UILabel("days:______________________");

        if (!HEADLESS) {
            win.AddCol(0, new UILabel("Cells"));
            win.AddCol(1, days);
            win.AddCol(0, Cell_vis);
            win.AddCol(2, new UILabel("CXCL9"));
            win.AddCol(2, CXCL9_vis);
            win.AddCol(2, new UILabel("TGF-\u03B2"));
            win.AddCol(2, TGFB_vis);
            win.RunGui();
        }

        // -------------------------
        // GIF makers
        // -------------------------
        GifMaker gm_Cell_vis = new GifMaker(subfolder.concat("/CellVid.gif"), 100, true);
        GifMaker gm_CXCL9_vis = new GifMaker(subfolder.concat("/CXCL9.gif"), 100, true);
        GifMaker gm_TGFB_vis = new GifMaker(subfolder.concat("/").concat("TGFBVid").concat(".gif"), 100, true);

        // -------------------------
        // Bone file selection and grid creation
        // -------------------------
        String Bone_Filename = null;
        if (LOCAL) {
            Bone_Filename = "/Users/80024703/Desktop/code/Bone/BAout_2020May5_Sim14.csv";
            //Bone_Filename = "/Users/80024703/Desktop/SmallBone.csv";
            //Bone_Filename = "/Users/80024703/Desktop/bone_data/day_269.csv";
            //Bone_Filename = "/Users/80024703/Desktop/bone_data/day_394.csv";
        } else {
            Bone_Filename = "Bone/BAout_2020May5_Sim14.csv";
        }

        // Create grid and RNG (give reproducible seed if you want deterministic runs)
        Rand rng = new Rand(simID + 1000L * prow); // deterministic seed per sim
        simpleBoneGrid g = new simpleBoneGrid(xDim, yDim, rng, Bone_Filename);


        // Set parameters from param_list if param sweep active
        if (PARAM_SWEEP && param_list != null) {
            g.SetParams(prow, param_list);
        }

        // Record Output (create PopOut.csv in this sim folder)
        g.newFileIO(subfolder, "w");

        // Initialize model
        g.InitBone( condition[0]); // First condition is myeloma cells

        // List for initializing vessels
        List<Integer> vessLocations = new ArrayList<>();

        double reducedMeanFraction = 0.2;
        int totalTcells = condition[1];
        int activeTcells = (int) Math.floor((1 - reducedMeanFraction) * totalTcells);
        int activeExtcells = totalTcells - activeTcells; // remainder to keep sum exact
        int cd8PopMax = activeTcells;
        int lowerCD8Thresh = (int) Math.floor(activeTcells * .46);

        // treg initial conditions
        double tregFraction = .11;
        double tregPercent = (((0.6 * tregFraction) / 100.0) * (xDim * yDim));
        int activeTregs = (int) Math.floor(tregPercent);
        int maxTregs = 30;

        // exhaustion parameters
        double dailyExhaustionProb = 0.00131;
        double dailyTcellDecreaseProb = 0.00333;
        double dailyTregIncreasenProb = 0.000545;

        // cd8 tcell initial conditions
        double tce_reducedMeanFraction = 0.2;
        int tce_totalTcells = condition[1];
        int tce_activeTcells = (int) Math.floor((1 - tce_reducedMeanFraction) * tce_totalTcells);
        int tce_activeExtcells = tce_totalTcells - tce_activeTcells; // remainder to keep sum exact
        int tce_cd8PopMax = tce_activeTcells;
        int tce_lowerCD8Thresh = (int) Math.floor(tce_activeTcells * .46);

        // treg initial conditions
        double tce_tregFraction = .11;
        double tce_tregPercent = (((0.6 * tce_tregFraction) / 100.0) * (xDim * yDim));
        int tce_activeTregs = (int) Math.floor(tce_tregPercent);
        int tce_maxTregs = 30;

        // exhaustion parameters
        double tce_dailyExhaustionProb = 0.00131;
        double tce_dailyTcellDecreaseProb = 0.00333;
        double tce_dailyTregIncreasenProb = 0.000545;


        // cd8 tcell initial conditions
        double no_tce_reducedMeanFraction = 0.20;
        int no_tce_totalTcells = (int) Math.floor(condition[1]*.20);
        int no_tce_activeTcells = (int) Math.floor((1 - no_tce_reducedMeanFraction) * no_tce_totalTcells);
        int no_tce_activeExtcells = no_tce_totalTcells - no_tce_activeTcells; // remainder to keep sum exact
        int no_tce_cd8PopMax = no_tce_activeTcells;
        int no_tce_lowerCD8Thresh = (int) Math.floor(no_tce_activeTcells * .46);

        // treg initial conditions
        double no_tce_tregFraction = .11;
        double no_tce_tregPercent = (((0.6 * no_tce_tregFraction) / 100.0) * (xDim * yDim));
        int no_tce_activeTregs = (int) Math.floor(no_tce_tregPercent);
        int no_tce_maxTregs = 10;

        // exhaustion parameters
        double no_tce_dailyExhaustionProb = 0.0;
        double no_tce_dailyTcellDecreaseProb = 0.0;
        double no_tce_dailyTregIncreasenProb = 0.0;

        boolean initial_recruitment = false;
        boolean initial_vessel = false;

        int timestepsPerDay = 240;
        int therapyPeriodWeeks = 1;
        // --- time loop (drop-in replacement) ---
        for (int i = 0; i < g.Nts; i++) {
            double[] Cell_Counts = g.CellCounts();
            timeStep = i;
            // Import external myeloma snapshots (unchanged)
            if (daysPassed > 0 & runPar & runPareShare) {
                int maxImport = 1; // adjust for how many per timestep
                for (int k = 0; k < maxImport; k++) {
                    CellSnapshot snap = myelomaTransferQueue.poll(); // thread-safe
                    if (snap != null) {
                        g.importSnapshotToGrid(snap);
                    } else {
                        break; // queue empty
                    }
                }
            }

            if (!HEADLESS) {
                win.TickPause(10);
            }


            if (Cell_Counts[7] + Cell_Counts[8] + Cell_Counts[9] <= 0 && i % 240 == 0 && CIRCULATING) {
                Random rand = new Random();
                int j = 0; // Counter for successfully placed T-cells
                int maxAttempts = 1000; // Prevent infinite loops
                int attempts = 0; // Counter for total attempts
                Collections.shuffle(vessLocations, rand);
                int[] boundaries = g.BoundaryIs();

                while (j < maxAttempts) {
                    attempts++;

                    // Try placing near a vessel
                    boolean placed = false;
                    int randomVessLocation = vessLocations.get(rand.nextInt(vessLocations.size()));

                    if (g.GetAgent(randomVessLocation).type == bloodVessel) {
                        int[] movdivHood = MooreHood(true);
                        int emptyNeighbors = g.MapEmptyHood(movdivHood, randomVessLocation);

                        if (emptyNeighbors > 0) {
                            int chosenIndex = g.rn.Int(emptyNeighbors);
                            int chosenCell = movdivHood[chosenIndex];

                            if (g.GetAgent(chosenCell) == null) {
                                simpleBoneCell c = g.NewAgentSQ(chosenCell);
                                c.type = MM;

                                j++;
                                placed = true;
                            }
                        }
                    }

                    // If placement near vessel failed, fall back to random boundary
                    if (!placed) {
                        int tries = 0;
                        while (tries < boundaries.length) {
                            int randomBoundaryIdx = rand.nextInt(boundaries.length);
                            int candidate = boundaries[randomBoundaryIdx];
                            if (g.GetAgent(candidate) == null) {
                                simpleBoneCell c = g.NewAgentSQ(candidate);
                                c.type = MM;

                                j++;
                                break;
                            }
                            tries++;
                        }
                    }

                    if (attempts >= maxAttempts) {
                        break;
                    }
                }
            }


            // compute days/weeks for scheduling (use same day-length logic as your daily recording)
            int daysPassedCalc = (int) Math.floor(i / (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)));
            int weeksPassed = daysPassedCalc / 7;

            // Weekly therapy toggle (mirrors your commented logic)
            boolean TCE_active_this_week = ((weeksPassed / therapyPeriodWeeks) % 2 == 0);

            // Switch TCELL at start of each week (only when we cross into a new day boundary)
            if (daysPassedCalc % 7 == 0 && i % (int)(24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0 && BIWEEKLY) {
                reducedMeanFraction = tce_reducedMeanFraction;
                totalTcells = tce_totalTcells;
                activeTcells = tce_activeTcells;
                activeExtcells = tce_activeExtcells; // remainder to keep sum exact


                // treg initial conditions
                tregFraction = tce_tregFraction;
                activeTregs = tce_activeTregs;
                maxTregs = tce_maxTregs;

                // exhaustion parameters
                dailyExhaustionProb = tce_dailyExhaustionProb;
                dailyTcellDecreaseProb = tce_dailyTcellDecreaseProb;
                dailyTregIncreasenProb = tce_dailyTregIncreasenProb;

            }
            else if (BIWEEKLY) {
                reducedMeanFraction = no_tce_reducedMeanFraction;
                totalTcells = no_tce_totalTcells;
                activeTcells = no_tce_activeTcells;
                activeExtcells = no_tce_activeExtcells; // remainder to keep sum exac

                // treg initial conditions
                tregFraction = no_tce_tregFraction;
                activeTregs = no_tce_activeTregs;
                maxTregs = no_tce_maxTregs;

                // exhaustion parameters
                dailyExhaustionProb = 0;
                dailyTcellDecreaseProb = 0;
                dailyTregIncreasenProb = 0;
            }


            if ((Cell_Counts[7] + Cell_Counts[8] + Cell_Counts[9]) < condition[0] * .25 && ADAPTIVE){
                reducedMeanFraction = tce_reducedMeanFraction;
                totalTcells = tce_totalTcells;
                activeTcells = tce_activeTcells;
                activeExtcells = tce_activeExtcells; // remainder to keep sum exact
                cd8PopMax = tce_cd8PopMax;
                lowerCD8Thresh = tce_lowerCD8Thresh;

                // treg initial conditions
                tregFraction = tce_tregFraction;
                activeTregs = tce_activeTregs;
                maxTregs = tce_maxTregs;

                // exhaustion parameters
                dailyExhaustionProb = tce_dailyExhaustionProb;
                dailyTcellDecreaseProb = tce_dailyTcellDecreaseProb;
                dailyTregIncreasenProb = tce_dailyTregIncreasenProb;
            }
            else if ((Cell_Counts[7] + Cell_Counts[8] + Cell_Counts[9] >= condition[0] * .25 && ADAPTIVE)){
                reducedMeanFraction = no_tce_reducedMeanFraction;
                totalTcells = no_tce_totalTcells;
                activeTcells = no_tce_activeTcells;
                activeExtcells = no_tce_activeExtcells; // remainder to keep sum exact
                cd8PopMax = no_tce_cd8PopMax;
                lowerCD8Thresh = no_tce_lowerCD8Thresh;

                // treg initial conditions
                tregFraction = no_tce_tregFraction;
                activeTregs = no_tce_activeTregs;
                maxTregs = no_tce_maxTregs;

                // exhaustion parameters
                dailyExhaustionProb = 0;
                dailyTcellDecreaseProb = 0;
                dailyTregIncreasenProb = 0;
            }

            if (TCELL && i == 0) {

                int InitTcells = 500;
                int j = 0;

                int maxAttempts = 100000;
                int attempts = 0;

                Random rand = new Random();
                int[] boundaries = g.BoundaryIs();

                while (j < InitTcells && attempts < maxAttempts) {

                    attempts++;

                    // -------- Try random location on grid --------
                    int chosenCell = rand.nextInt(g.length);

                    if (g.GetAgent(chosenCell) == null) {

                        simpleBoneCell c = g.NewAgentSQ(chosenCell);

                        if (g.rn.Double() < reducedMeanFraction) {
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(10, 1, 10, 34);
                        } else if (g.rn.Double() < .1) {
                            c.type = supressorTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        } else {
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(20, 1, 1, 22);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }

                        j++;
                        continue;
                    }



                    // -------- Fallback to boundary --------
                    int randomBoundaryIdx = rand.nextInt(boundaries.length);
                    chosenCell = boundaries[randomBoundaryIdx];

                    if (g.GetAgent(chosenCell) == null) {

                        simpleBoneCell c = g.NewAgentSQ(chosenCell);

                        if (g.rn.Double() < reducedMeanFraction) {
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(10, 1, 10, 34);
                        } else if (g.rn.Double() < .1) {
                            c.type = supressorTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        } else {
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(20, 1, 1, 22);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }
                        j++;
                    }
                }
            }

            if (TCELL && i < 0) {

                int InitTcells = 2;
                int j = 0;

                int maxAttempts = 100000;
                int attempts = 0;

                Random rand = new Random();
                int[] boundaries = g.BoundaryIs();

                while (j < InitTcells && attempts < maxAttempts) {

                    attempts++;

                    // -------- Try random location on grid --------
                    int chosenCell = rand.nextInt(g.length);

                    if (g.GetAgent(chosenCell) == null) {

                        simpleBoneCell c = g.NewAgentSQ(chosenCell);

                        if (g.rn.Double() < reducedMeanFraction) {
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(10, 1, 10, 34);
                        } else {
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(20, 1, 1, 22);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }

                        j++;
                        continue;
                    }



                    // -------- Fallback to boundary --------
                    int randomBoundaryIdx = rand.nextInt(boundaries.length);
                    chosenCell = boundaries[randomBoundaryIdx];

                    if (g.GetAgent(chosenCell) == null) {

                        simpleBoneCell c = g.NewAgentSQ(chosenCell);

                        if (g.rn.Double() < reducedMeanFraction) {
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(10, 1, 10, 34);
                        } else {
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(20, 1, 1, 22);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }
                        j++;
                    }
                }
            }

            // Model step and drawing (same as before)
            g.ModelStep(i, Cell_Counts, simID);

            if (!HEADLESS) {
                g.Draw(Cell_vis, days, i, simID);
                g.DrawCXCL9(CXCL9_vis);
                g.DrawTGFB(TGFB_vis);
            }
            // daily recording: your original condition:
            if (i % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                g.RecordOut(g.out, i, TREATMENT_ON, BORTEZOMIB, MYELOMA);
                g.RecordClones(g.clones, i);
                g.RecordBones(g.bones, i);
                g.RecordLocs(g.locations, i);
                g.RecordTGFB(g.TGFBout, i);
                g.RecordParamsOut(g.paramsOut,xDim, yDim, reducedMeanFraction, totalTcells, activeTcells, activeExtcells,
                        tregFraction, activeTregs, maxTregs, dailyExhaustionProb, dailyTcellDecreaseProb,
                        dailyTregIncreasenProb, timestepsPerDay, therapyPeriodWeeks);

                daysPassed += 1;

                // Export myeloma cells for transfer
                for (simpleBoneCell c : g) {
                    if (c.type == MM & runPareShare) {
                        if (g.rn.Double() < 0.01) { // tune probability
                            CellSnapshot snap = g.exportSnapshot(c);
                            myelomaTransferQueue.add(snap);
                            c.Dispose(); // remove from original grid
                        }
                    }
                }
            }
        } // end time loop


        // -------------------------
        // cleanup: close files, gifs, window
        // -------------------------
        g.closeFileIO();
        gm_Cell_vis.Close();
        gm_CXCL9_vis.Close();
        gm_TGFB_vis.Close();

        if (!HEADLESS) {
            win.Close();
        }
    }


    //sample from a bounded  distribution
    public double boundedGaussian(double mean, double dev, double min, double max) {
        double gauss = rn.Gaussian(0, 1);
        double val = dev * gauss + mean;
        while (val > max || val < min) {
            gauss = rn.Gaussian(0, 1);
            val = dev * gauss + mean;
        }
        return val;
    }

    public void newFileIO (String projPath, String mode) {

        out = new FileIO(projPath + "PopOut.csv", mode);
//        params = new FileIO(projPath + "params.csv",mode);
        clones = new FileIO(projPath + "clones.csv", mode);

        bones = new FileIO(projPath + "boneOut.csv", mode);

        locations = new FileIO(projPath + "cellLocs.csv", mode);

        TGFBout=new FileIO(projPath + "TGFBOut.csv", mode);

        paramsOut = new FileIO(projPath + "paramsOut.txt",mode);


        if(mode=="w") {

            out.Write("Timestep" + "," + "BONE" + "," + "pOB" + "," + "aOB" + "," + "pOC" + "," + "aOC" + "," + "MSC" + "," + "LINING" + "," + "S_MM" + "," + "R_MM"+"," +"AL_MM"+ ","+"TCell"+","+"ExtTcell"+"," +"T-reg"+"," + "Naive Tcell" +","+ "TREATMENT_ON" + "," + "BORTEZOMIB" + "," + "MYELOMA" + "\n");
            clones.Write("Timestep" + "," + "SimID" + "," + "cloneID" + ","+"cellID"+","+ "MHCI" + "," + "BCMA" + "\n");
            bones.Write("Timestep" + "," + "SimID" + "," + "Position" + "\n");
            locations.Write("Timestep" + "," + "SimID" + "," + "Position" + ","+"Type"+"\n");
            TGFBout.Write("Timestep"+","+ "TGFB" + "\n");
            paramsOut.Write("Header");
        }

    }

    public void closeFileIO () {

        out.Close();
        clones.Close();
        bones.Close();
        locations.Close();
        TGFBout.Close();
        paramsOut.Close();
    }

    public void SetParams(int prow, ArrayList<String> param_list){
        //returns an array list of all lines from the file as stringsftype == lining

        String[] split_param_list = param_list.get(prow).split(",");


        MYELOMA = Boolean.parseBoolean(split_param_list[0]);
        BORTEZOMIB = Boolean.parseBoolean(split_param_list[1]);
        pmutate = Double.parseDouble(split_param_list[2]);
        EMDR = Boolean.parseBoolean(split_param_list[3]);
    }

    // Helper method to check if a cell is in the bone lining area
    private boolean isInLining(int x, int y) {
        simpleBoneCell cell = GetAgent(x, y);
        return (cell != null && cell.type == LINING);
    }
    public void InitBone(int initMyeloma) {


//  FOR IRREGULAR BONE
        int xinit, yinit;
        ArrayList<String> input_data = InitialBone.Read();
        String[] split_input_data =input_data.get(0).split(",");

        //Place bone
        for (int index=1; index<split_input_data.length; index++){
            NewAgentSQ(Integer.parseInt(split_input_data[index])).type=BONE;
            GetAgent(Integer.parseInt(split_input_data[index])).Init();
            InitBoneList.add(Integer.parseInt(split_input_data[index]));
            AllBoneList.add(GetAgent(Integer.parseInt(split_input_data[index])));
        }
        for (int index=1; index<split_input_data.length; index++){
            if(GetAgent(Integer.parseInt(split_input_data[index])).MarrowInHood()==true){
                GetAgent(Integer.parseInt(split_input_data[index])).type=LINING;
                GetAgent(Integer.parseInt(split_input_data[index])).liningAge = TURNOVER_TIME;
                LiningList.add(GetAgent(Integer.parseInt(split_input_data[index])));
            }
        }
        init_BA=InitBoneList.size();
        MarrowArea = (xDim*yDim)-init_BA;//(xDimBone*yDimBone); //0.12 Bone, 0.88 Marrow
        int myelomaCellsToPlace = initMyeloma;
        int placedMyelomaCells = 0;
        int boneProximityDistance = 10; // Maximum initial distance from bone
        double bcmaNegfraction = 0.0;
        int nextCloneID = 1;


        int attempts = 0;
        int maxAttempts = 100000; // safety to prevent infinite loop

        while (placedMyelomaCells < myelomaCellsToPlace && attempts < maxAttempts) {

            int x = rn.Int(xDim);
            int y = rn.Int(yDim);

            // Skip if occupied
            if (PopAt(x, y) != 0) {
                attempts++;
                continue;
            }

            // Check if near bone
            boolean isNearBone = false;

            for (int xi = Math.max(0, x - boneProximityDistance);
                 xi <= Math.min(xDim - 1, x + boneProximityDistance); xi++) {

                for (int yi = Math.max(0, y - boneProximityDistance);
                     yi <= Math.min(yDim - 1, y + boneProximityDistance); yi++) {

                    if (GetAgent(xi, yi) != null && GetAgent(xi, yi).type == BONE) {

                        double distance = Math.sqrt(
                                Math.pow(xi - x, 2) +
                                        Math.pow(yi - y, 2));

                        if (distance <= boneProximityDistance) {
                            isNearBone = true;
                            break;
                        }
                    }
                }
                if (isNearBone) break;
            }

            if (!isNearBone) {
                attempts++;
                continue;
            }

            // Place myeloma cell
            NewAgentSQ(x, y).type = MM;
            GetAgent(x, y).bcmaExpression = 1;
            GetAgent(x, y).mhcIExpression = 1;
            GetAgent(x, y).cloneID = nextCloneID++;
            GetAgent(x, y).cellID = GetAgent(x, y).cloneID;

            if (rn.Double() < bcmaNegfraction) {
                GetAgent(x, y).bcmaLoss = true;
                GetAgent(x, y).bcmaExpression = rn.Double();
            }

            placedMyelomaCells++;
        }




    }

    public void ModelStep(int time, double [] Cell_Counts, int simID) {

        //STEP 0: UPDATE GRIDTICK
        GridTick();
        int i=0;
        double stol = 1.0e-6;//1e-6; //steady-state tolerance

        //TGFB Diffusion coefficient
        Grid2Ddouble TGFB_xDiffArray = new Grid2Ddouble(xDim,yDim);
        Grid2Ddouble TGFB_yDiffArray = new Grid2Ddouble(xDim,yDim);

        TGFB_xDiffArray.SetAll(TGFB_DiffCoef);
        TGFB_yDiffArray.SetAll(TGFB_DiffCoef);

        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    TGFB_xDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (x!=xDim-1 && GetAgent(x+1, y) != null && (GetAgent(x+1, y).type == BONE || GetAgent(x+1, y).type == LINING)){
                    TGFB_xDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    TGFB_yDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (y!=yDim-1 && GetAgent(x, y+1) != null && (GetAgent(x, y+1).type == BONE || GetAgent(x, y+1).type == LINING)){
                    TGFB_yDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
            }
        }


        for (int x = 0; x < CXCL9.xDim; x++) {
            for (int y = 0; y < CXCL9.yDim; y++) {
                if (GetAgent(x,y)!=null && GetAgent(x, y).type == MM && GetAgent(x, y).bcmaExpression >0)  {
                    CXCL9.Add(x, y, CXCL9_productionRate/maxCXCL9 );
                } if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    CXCL9.Set(x, y, 0.1 * CXCL9_DiffCoef); //Originally had CXCL9_DiffCoef*0.1
                } else if (x!=xDim-1 && GetAgent(x+1, y) != null && (GetAgent(x+1, y).type == BONE || GetAgent(x+1, y).type == LINING)){
                    CXCL9.Set(x, y, 0.1 * CXCL9_DiffCoef); //Originally had CXCL9_DiffCoef*0.1
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    CXCL9.Set(x, y, 0.1 * CXCL9_DiffCoef); //Originally had CXCL9_DiffCoef*0.1
                } else if (y!=yDim-1 && GetAgent(x, y+1) != null && (GetAgent(x, y+1).type == BONE || GetAgent(x, y+1).type == LINING)){
                    CXCL9.Set(x, y, 0.1 * CXCL9_DiffCoef); //Originally had CXCL9_DiffCoef*0.1
                }
            }
        }

        //CXCL9 Diffusion
        CXCL9.DiffusionADI(CXCL9_DiffCoef);

        //Natural Decay of CXCL9
        CXCL9.MulAll(CXCL9_decayRate);
        CXCL9.Update();

        do {
            IterateTGFB(TGFB_xDiffArray,TGFB_yDiffArray);
            //Iterate
            i++;
        } while (i<N_TIMESTEP_PDE && MaxTdiff >= stol);

        if(tmax<TGFB.GetMax()){
            tmax=TGFB.GetMax();
        }

        /////////////////////////////////
        //STEP 3: ITERATE THROUGH CELLS//
        /////////////////////////////////

        CleanShuffle(rn);
        //ShuffleAgents(rn);
        for (simpleBoneCell c: this) {
            c.CellStep(time, Cell_Counts, simID);
        }


    }

    ////////////////////////////
    //Full Domain (No zoom-in)//
    ////////////////////////////
    public void Draw(UIGrid vis, UILabel days, int i, int simID) {
        days.SetText("days: "+i*convert_to_days);
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x,y);
                if (drawMe != null && drawMe.type==MM && drawMe.RESISTANT==true){
                    vis.SetPix(x,y, BLACK);
                } else if (drawMe != null && drawMe.type==MM && drawMe.bcmaLoss==true){
                    vis.SetPix(x,y, BLACK);
                }else if (drawMe != null && drawMe.type==MM && drawMe.simulationID!=simID) {
                    vis.SetPix(x, y, RGB(255, 255, 0));
                }else if (drawMe != null) {
                    //vis.SetPix(x, y, drawMe.color);
                    vis.SetPix(x, y, drawMe.type);
                } else{
                    vis.SetPix(x,y, RGB256(240, 220, 220)); //MARROW=LIGHT PINK
                }
            }
        }
        vis.SetString("Day: "+(int)(i*convert_to_days),1,yDim-1,BLACK, RGB256(240, 220, 220));

    }
    public void DrawTGFB(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else
                    vis.SetPix(x, y, HeatMapRGB(TGFB.Get(x, y)));

            }
        }
    }
    public void DrawCXCL9(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));

            }
        }
    }
    public void DrawIFNG(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else
                    vis.SetPix(x, y, HeatMapRGB(IFNG.Get(x, y)));

            }
        }
    }

    public void RecordOut(FileIO writeHere,int time, boolean treatment_on, boolean btz, boolean myeloma){
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int[] cts = new int[14];

        for (simpleBoneCell c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            } else if(c.type==pOB){
                //ct_pOB++;
                cts[1]++;
            } else if(c.type==aOB){
                //ct_aOB++;
                cts[2]++;
            } else if(c.type==pOC){
                //ct_pOC++;
                cts[3]++;
            } else if(c.type==aOC){
                //ct_aOC++;
                cts[4]++;
            } else if(c.type==MSC){
                //ct_MSC++;
                cts[5]++;
            } else if(c.type==LINING){
                //ct_LINING++;
                cts[6]++;
            } else if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
                cts[7]++;
            } else if(c.type==MM && c.RESISTANT) {
                cts[8]++;
            } else if(c.type==MM && c.bcmaLoss) {
                cts[9]++;
            } else if(c.type==activeTcell) {
                cts[10]++;
            } else if(c.type==EXHT_CELL) {
                cts[11]++;
            }
            else if(c.type == supressorTcell){
                cts[12]++;
            }
            else if (c.type == naiveTcell){
                cts[13]++;
            }

        }
        //population of one timestep per line
        writeHere.Write(time+",");
        writeHere.WriteDelimit(cts,",");
        writeHere.Write("," + treatment_on + "," + btz + "," + myeloma + "\n");
    }
    public void RecordClones(FileIO writeHere,int time){

        for (simpleBoneCell c : this) {
            if(c.type==MM) {
                int simID = c.simulationID;
                int cloneID = c.cloneID;
                int cellID = c.cellID;
                double mhcI = c.mhcIExpression;
                double bcma = c.bcmaExpression;
                writeHere.Write(time + "," + simID + "," + cloneID + ","+cellID+","+ mhcI + ","+ bcma + "\n");
            }
        }
    }

    public void RecordTGFB(FileIO writeHere, int time) {
        double sum = 0.0;
        int count = 0;

        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                sum += TGFB.Get(x, y);
                count++;
            }
        }

        double avgTGFB = sum / count;
        writeHere.Write(time + "," + avgTGFB + "\n");
    }
    public void RecordBones(FileIO writeHere,int time){

        for (simpleBoneCell c : this) {
            if(c.type==BONE) {
                int simID = c.simulationID;
                int position = c.Isq();
                writeHere.Write(time + "," + simID + "," + position + "\n");
            }
        }
    }

    public void RecordLocs(FileIO writeHere,int time){

        for (simpleBoneCell c : this) {
            int simID = c.simulationID;
            int position = c.Isq();
            String cellType = "UNKNOWN";
            if(c.type==pOB){
                cellType = "pOB";
            } else if(c.type==aOB){
                cellType = "aOB";
            } else if(c.type==pOC){
                cellType = "pOC";
            } else if(c.type==aOC){
                cellType = "aOC";
            } else if(c.type==MSC){
                cellType = "MSC";
            } else if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
                cellType = "Sensitive MM";
            } else if(c.type==MM && c.RESISTANT) {
                cellType = "Resistant MM";
            } else if(c.type==MM && c.bcmaLoss) {
                cellType = "AL_MM";
            } else if(c.type==activeTcell) {
                cellType = "Tcell";
            } else if(c.type==EXHT_CELL) {
                cellType = "Ext_Tcell";
            } else if(c.type == supressorTcell){
                cellType = "Treg";
            } else if (c.type == naiveTcell){
                cellType = "Naive Tcell";
            } else if(c.type == BONE){
                cellType = "BONE";
            }else if(c.type == LINING){
                cellType = "LINING";
            }
            writeHere.Write(time + "," + simID + "," + position +  "," + cellType +"\n");
        }
    }

    public void RecordParamsOut(FileIO writeHere,double xDim, double yDim, double reducedMeanFraction, double totalTcells, double activeTcells,
                                double activeExtcells, double tregFraction, double activeTregs, double maxTregs,
                                double dailyExhaustionProb, double dailyTcellDecreaseProb, double dailyTregIncreasenProb,
                                double timestepsPerDay, double therapyPeriodWeeks){
        // ---- WRITE ALL PARAMETERS TO paramsOut.txt ----
        writeHere.Write("=== MODEL PARAMETERS ===\n");

        writeHere.Write("xDim: " + xDim + "\n");
        writeHere.Write("yDim: " + yDim + "\n");
        writeHere.Write("reducedMeanFraction: " + reducedMeanFraction + "\n");
        writeHere.Write("totalTcells: " + totalTcells + "\n");
        writeHere.Write("activeTcells: " + activeTcells + "\n");
        writeHere.Write("activeExtcells: " + activeExtcells + "\n");
        writeHere.Write("tregFraction: " + tregFraction + "\n");
        writeHere.Write("activeTregs: " + activeTregs + "\n");
        writeHere.Write("maxTregs: " + maxTregs + "\n");

        writeHere.Write("dailyExhaustionProb: " + dailyExhaustionProb + "\n");
        writeHere.Write("dailyTcellDecreaseProb: " + dailyTcellDecreaseProb + "\n");
        writeHere.Write("dailyTregIncreasenProb: " + dailyTregIncreasenProb + "\n");

        writeHere.Write("timestepsPerDay: " + timestepsPerDay + "\n");
        writeHere.Write("therapyPeriodWeeks: " + therapyPeriodWeeks + "\n");

    }


    public double[] CellCounts(){
        double[] cts = new double[14];

        for (simpleBoneCell c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            } else if(c.type==pOB){
                //ct_pOB++;
                cts[1]++;
            } else if(c.type==aOB){
                //ct_aOB++;
                cts[2]++;
            } else if(c.type==pOC){
                //ct_pOC++;
                cts[3]++;
            } else if(c.type==aOC){
                //ct_aOC++;
                cts[4]++;
            } else if(c.type==MSC){
                //ct_MSC++;
                cts[5]++;
            } else if(c.type==LINING){
                //ct_LINING++;
                cts[6]++;
            } else if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
                cts[7]++;
            } else if(c.type==MM && c.RESISTANT){
                cts[8]++;
            } else if(c.type==MM && c.bcmaLoss) {
                cts[9]++;
            } else if(c.type == activeTcell){
                cts[10]++;
            } else if(c.type == EXHT_CELL){
                cts[11]++;
            } else if(c.type == supressorTcell){
                cts[12]++;
            }
            else if (c.type == naiveTcell){
                cts[13]++;
            }
        }
        return cts;
    }



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      MAIN                                                      //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static void main(String[] args) throws InterruptedException {
        if (HEADLESS) {
            System.setProperty("java.awt.headless", "true");
        }

        // base folder and timestamp (this is the code you kept before)
        String sdf = new SimpleDateFormat("yyyyMMMd").format(new Date());
        String fn = "Bone_" + sdf;
        File baseDir = new File(fn);
        baseDir.mkdir();

        // PARAM SWEEP reading
        int param_list_size;
        ArrayList<String> param_list = null;
        if (PARAM_SWEEP) {
            FileIO Params = new FileIO("Bone/boneRemodeling_2022May17/params.csv", "r");
            param_list = Params.Read();
            param_list_size = param_list.size();
        } else {
            param_list_size = 2;
        }

        // Thread pool (limit threads if you want)
        int numThreads = Math.min(numSims, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        int[][] initialConditions = {
                {3000, 200}
        };
        if (runPar){
            for (int prow = 1; prow < param_list_size; prow++) {
                for (int sim = 0; sim < numSims; sim++) {
                    final int simID = sim;
                    final int row = prow;
                    final String projectBase = fn;
                    final ArrayList<String> paramsCopy = param_list; // may be null

                    final int[] condition = initialConditions[simID % initialConditions.length];

                    executor.submit(() -> {
                        runSimulation(simID, row, projectBase, paramsCopy, condition);
                    });
                }
            }
            executor.shutdown();
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            System.out.println("All simulations finished.");
        }
        else{
            for (int prow = 1; prow < param_list_size; prow++) {
                for (int sim = 0; sim < numSims; sim++) {
                    final int simID = sim;
                    final int row = prow;
                    final String projectBase = fn;
                    final ArrayList<String> paramsCopy = param_list; // may be null
                    final int[] condition = initialConditions[simID % initialConditions.length];
                    runSimulation(simID, row, projectBase, paramsCopy, condition);
                }
            }
        }


    }



}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     CELL CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class simpleBoneCell extends AgentSQ2Dunstackable<simpleBoneGrid> {

    ///////////////
    //CELL FIELDS//
    ///////////////

    public int type;

    public int color;

    //Lining
    int liningAge = 0; //counter for when lining old enough for remodeling event
    //MM
    boolean RESISTANT = false;

    boolean bcmaLoss = false;
    boolean TGFB_on = false;
    double tcellAge = 0;
    double lifespan = 0;

    double pd_1 = 0;
    double pd_l1 = 0;
    boolean myeloma_bound = false;
    boolean mhcLoss = false;
    double bcmaExpression = 1.0;
    double mhcIExpression = 1.0;
    int simulationID;
    int cloneID;
    double p_kill = ProbScale(.25, TIMESTEP_AGENT);
    double extp_kill = ProbScale(.1, TIMESTEP_AGENT);

    int cellID;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CELL METHODS//////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. CellStep:                                                                                        //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    public int seekCXCL9() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] CXCL9_levels = new double[9]; // stores CXCL9 levels at the 8 directions and the center

        // Extracting CXCL9 levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.CXCL9.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); // cumulative probabilities

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]); // add index to list

                double P = 0;
                switch (i) {
                    case 1: // right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 2: // left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 3: // up
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 4: // down
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 5: // top-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 6: // top-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 7: // bottom-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                    case 8: // bottom-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                }

                P = Math.max(P, 0); // ensure non-negative probability

                if (cProbArray.size() > 0) {
                    cProbArray.add(P + cProbArray.get(cProbArray.size() - 1));
                } else {
                    cProbArray.add(P);
                }
                ProbSum += P;
            }
        }

        int moveToIndex = G.tmoveHood[G.rn.Int(neighbors)]; // Default to the center index if no movement

        for (int i = 0; i < cProbArray.size(); i++) {
            if (rnum <= cProbArray.get(i) / ProbSum) {
                moveToIndex = emptyHood.get(i);
                break;
            }
        }

        return moveToIndex;
    }

    public int dysseekCXCL9() {
        double exhaustionLevel = 1;

        int neighbors = MapHood(G.tmoveHood);
        double[] CXCL9_levels = new double[9];

        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.CXCL9.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;

        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>();

        // ↓↓↓ Exhaustion modifies behavior ↓↓↓

        double effectiveTaxis = G.Tcell_TaxisCoeff * (1.0 - exhaustionLevel);

        // Noise increases with exhaustion
        double noiseStrength = exhaustionLevel * G.Tcell_TaxisCoeff;

        // Probability of fully random step
        double randomMoveProb = 0.2 * exhaustionLevel;

        // If exhausted cell ignores gradient entirely
        if (G.rn.Double() < randomMoveProb) {
            return G.tmoveHood[G.rn.Int(neighbors)];
        }

        for (int i = 0; i < neighbors; i++) {

            if (G.GetAgent(G.tmoveHood[i]) == null) {

                emptyHood.add(G.tmoveHood[i]);

                double P = 0;
                double gradientTerm = 0;

                switch (i) {
                    case 1:
                        gradientTerm = (CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 2:
                        gradientTerm = -(CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 3:
                        gradientTerm = (CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 4:
                        gradientTerm = -(CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 5:
                        gradientTerm = (CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 6:
                        gradientTerm = -(CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 7:
                        gradientTerm = (CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                    case 8:
                        gradientTerm = -(CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                }

                // Add exhaustion-dependent noise
                double noisyGradient = gradientTerm +
                        G.rn.Double() * noiseStrength;

                P = G.Tcell_DiffCoef +
                        (effectiveTaxis * G.maxCXCL9 / 8.0) * noisyGradient;

                P = Math.max(P, 0);

                if (cProbArray.size() > 0) {
                    cProbArray.add(P + cProbArray.get(cProbArray.size() - 1));
                } else {
                    cProbArray.add(P);
                }

                ProbSum += P;
            }
        }

        int moveToIndex = G.tmoveHood[G.rn.Int(neighbors)];

        if (ProbSum > 0) {
            for (int i = 0; i < cProbArray.size(); i++) {
                if (rnum <= cProbArray.get(i) / ProbSum) {
                    moveToIndex = emptyHood.get(i);
                    break;
                }
            }
        }

        return moveToIndex;
    }

    public int seekTGFBDIAG() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] IFNG_levels = new double[9]; // stores IFNG levels at the 8 directions and the center

        // Extracting IFNG levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            IFNG_levels[i] = G.TGFB.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); // cumulative probabilities

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]); // add index to list

                double P = 0;
                switch (i) {
                    case 1: // right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[1] - IFNG_levels[2]);
                        break;
                    case 2: // left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[1] - IFNG_levels[2]);
                        break;
                    case 3: // up
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[3] - IFNG_levels[4]);
                        break;
                    case 4: // down
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[3] - IFNG_levels[4]);
                        break;
                    case 5: // top-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[5] - IFNG_levels[6]);
                        break;
                    case 6: // top-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[5] - IFNG_levels[6]);
                        break;
                    case 7: // bottom-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[7] - IFNG_levels[8]);
                        break;
                    case 8: // bottom-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (IFNG_levels[7] - IFNG_levels[8]);
                        break;
                }

                P = Math.max(P, 0); // ensure non-negative probability

                if (cProbArray.size() > 0) {
                    cProbArray.add(P + cProbArray.get(cProbArray.size() - 1));
                } else {
                    cProbArray.add(P);
                }
                ProbSum += P;
            }
        }

        int moveToIndex = G.tmoveHood[G.rn.Int(neighbors)]; // Default to the center index if no movement

        for (int i = 0; i < cProbArray.size(); i++) {
            if (rnum <= cProbArray.get(i) / ProbSum) {
                moveToIndex = emptyHood.get(i);
                break;
            }
        }

        return moveToIndex;
    }

    public int seekTGFB(double DiffCoef, double TaxisCoef) {
        int neighbors = MapHood(G.moveHood); //includes self
        double TGFB_this = G.TGFB.Get(G.moveHood[0]); //or G.TGFB.Get(this.Xsq(), this.Ysq());
        double TGFB_right = G.TGFB.Get(G.moveHood[1]); //or G.TGFB.Get(this.Xsq() + 1, this.Ysq());
        double TGFB_left = G.TGFB.Get(G.moveHood[2]); //or G.TGFB.Get(this.Xsq() - 1, this.Ysq());
        double TGFB_up = G.TGFB.Get(G.moveHood[3]); //or G.TGFB.Get(this.Xsq(), this.Ysq() + 1);
        double TGFB_down = G.TGFB.Get(G.moveHood[4]);// or G.TGFB.Get(this.Xsq(), this.Ysq() - 1);
        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); //cumulative probabilities


        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.moveHood[i]) == null) {
                emptyHood.add(G.moveHood[i]); //add index to list

                switch (i) {
                    case 0: //no movement
                        double P0 = 1 - 4 * DiffCoef - TaxisCoef * maxTGFB * (TGFB_right + TGFB_left - 4 * TGFB_this + TGFB_up + TGFB_down);
                        if (P0 < 0) {
                            P0 = 0;
                            //System.out.println("Warning: P0<0");
                        }
                        cProbArray.add(P0);
                        ProbSum += P0;
                        break;
                    case 1: //right
                        double P2 = DiffCoef + (TaxisCoef * maxTGFB) / 4 * (TGFB_right - TGFB_left);
                        if (P2 < 0) {
                            P2 = 0;
                            //System.out.println("Warning: P2<0");
                        }
                        if (cProbArray.size() > 0) {
                            cProbArray.add(P2 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P2);
                        }
                        ProbSum += P2;
                        break;
                    case 2: //left
                        double P1 = DiffCoef - (TaxisCoef * maxTGFB) / 4 * (TGFB_right - TGFB_left);
                        if (P1 < 0) {
                            P1 = 0;
                            //System.out.println("Warning: P1<0");
                        }
                        if (cProbArray.size() > 0) {
                            cProbArray.add(P1 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P1);
                        }
                        ProbSum += P1;
                        break;
                    case 3: //up
                        double P4 = DiffCoef + (TaxisCoef * maxTGFB) / 4 * (TGFB_up - TGFB_down);
                        if (P4 < 0) {
                            P4 = 0;
                            //System.out.println("Warning: P4<0");
                        }
                        if (cProbArray.size() > 0) {
                            cProbArray.add(P4 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P4);
                        }
                        ProbSum += P4;
                        break;
                    case 4: //down
                        double P3 = DiffCoef - (TaxisCoef * maxTGFB) / 4 * (TGFB_up - TGFB_down);
                        if (P3 < 0) {
                            P3 = 0;
                            //System.out.println("Warning: P3<0");
                        }
                        if (cProbArray.size() > 0) {
                            cProbArray.add(P3 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P3);
                        }
                        ProbSum += P3;
                        break;
                }

            }
        }

        // If no valid movement probabilities, stay put
        if (ProbSum <= 0 || cProbArray.size() == 0) {
            return G.moveHood[0];
        }

        // Scale random number to total probability
        double scaledRand = G.rn.Double() * ProbSum;

        // Select move based on cumulative distribution
        for (int j = 0; j < cProbArray.size(); j++) {
            if (scaledRand <= cProbArray.get(j)) {
                return emptyHood.get(j);
            }
        }

        // Fallback (should not hit)
        return G.moveHood[0];

    }

    void Init() {
        if (type == BONE || type == LINING) { //Include aOB because they will turn into bone
            G.count_BA++;
        }
    }

    public double Prob_Divide(double TGFBval, double MaxDivRate, double halfmax) {
        double n = 2;
        return MaxDivRate * (1 / (1 + Math.pow((halfmax / TGFBval), n)));
    }

    public boolean MarrowInHood() {
        int[] MarrowHood = VonNeumannHood(false); //4 neighbors
        int options = MapEmptyHood(MarrowHood); //occupied positions surrounding pOB
        return options > 0;
    }

    void Tcell_Kill() {
        if (type == MM) {
            this.Dispose();
        }
    }

    public void CellStep(int time, double[] Cell_Counts, int simID) {


        ///////////
        //MYELOMA//
        ///////////
        if (type == MM) {

            if (time == 0) {

                this.simulationID = simID;
            }
//            if (this.simulationID!=simID) {
//                color = RGB(255, 255, 0);
//            }
            double rn_BirthDeath = G.rn.Double();
            double pdiv;
            double pdeath = 0;

            if (G.TGFB.Get(this.Isq()) < G.TGFBthresh) { // && rn_BirthDeath < G.MM_DEATH) {
                pdeath = G.MM_DEATH; //No BDF survival advantage (no tx)

            } else if (G.TGFB.Get(this.Isq()) >= G.TGFBthresh) { // && rn_BirthDeath < G.MM_DEATH_BDF){
                pdeath = G.MM_DEATH_BDF; //BDF survival advantage (no tx)

            }
            double scaleFactor = 0.7 + (0.3 * this.bcmaExpression); // bcma=0 → 0.7, bcma=1 → 1.0
            pdiv = Prob_Divide(G.TGFB.Get(this.Isq()), 1.0 / 1440 * MinToHour, Math.sqrt(3.0) * G.TGFBthresh); //same as pOB

            ///////////
            //MM dies//
            ///////////
            if (rn_BirthDeath < ProbScale(pdeath, TIMESTEP_AGENT)) {//rn_BirthDeath < pdeath){
                //MM dies
                color = WHITE;
                Dispose();
            }

            int[] Hood = MooreHood(true); // For division and movement
            int options = MapOccupiedHood(Hood); // Mapping occupied spots
            double x = 5.5e-4;
            for (int j = 0; j < options; j++) {
                if (G.GetAgent(Hood[j]) != null && G.GetAgent(Hood[j]).type == BONE || G.GetAgent(Hood[j]).type == LINING) {
                    if (G.rn.Double() < x) {
                        G.TGFBtimer[this.Xsq()][this.Ysq()] = 0;
                        TGFB_on = true;
                        G.GetAgent(Hood[j]).Dispose();
                    }
                    }else {
                    TGFB_on = false;
                }
            }
            if (rn_BirthDeath < ProbScale(pdiv, TIMESTEP_AGENT)) {

                //////////////
                //MM divides//
                //////////////
                // Preallocate a neighborhood array (outside loop, reuse for speed)
                int[] hood = MooreHood(true);
                int emptyNeighbors = MapEmptyHood(hood);
                if (emptyNeighbors > 0) {
                    // Choose a random empty spot from the mapped hood
                    int rIndex = G.rn.Int(emptyNeighbors);
                    int chosenCell = hood[rIndex];

                    if (G.GetAgent(chosenCell) == null) { // Confirm still empty
                        simpleBoneCell child = G.NewAgentSQ(chosenCell);
                        child.type = this.type;
                        child.RESISTANT = this.RESISTANT;
                        child.simulationID = this.simulationID;
                        child.cloneID = this.cloneID;
                        child.cellID = cellID++;

                        if (G.rn.Double() < G.antigenLoss && mhcLoss == false) {
                            child.mhcIExpression = 0;
                            child.cloneID = this.cloneID;
                            child.mhcLoss = true;
                        }

                        // Antigen loss mutation
                        if (G.rn.Double() < G.antigenLoss && bcmaLoss == false) {
                            // Pick one of four discrete expression levels: 0, 0.25, 0.5, or 0.75
                            //double[] possibleValues = {0.0, 0.25, 0.5, 0.75};
                            //int index = G.rn.Int(possibleValues.length); // random integer 0–3
                            //double newExpression = possibleValues[index];
                            if (child != null) {
                                child.bcmaExpression = 0;
                                child.bcmaLoss = true;
                            }
                        } else {
                            child.bcmaExpression = this.bcmaExpression;
                            child.mhcLoss = mhcLoss;
                        }
                    }
                }
            }
//            int[] hood = VonNeumannHood(true);
//            int emptyNeighbors = MapEmptyHood(hood);
//            if (emptyNeighbors > 0){
//            MoveSQ(seekTGFB(0,5.0e9*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP)));//only towards resorption
//                }

        }
        if (type == naiveTcell) {

            if(G.TGFB.Get(Isq()) >= G.TGFBthresh){
                this.pd_1+=1;
            }
            this.myeloma_bound = false;
            int[] hood = MooreHood(true);
            // Increment T cell age every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }
            // Age-related death logic
            if (this.tcellAge >= this.lifespan) {
                int emptyNeighbors = MapEmptyHood(hood);

                if (emptyNeighbors > 0 &&
                        G.rn.Double() < ProbScale(1/1440 *MinToHour, TIMESTEP_AGENT)) {

                    int chosenCell = hood[G.rn.Int(emptyNeighbors)];

                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                    child.pd_l1 = child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
                }
            }
            int emptyNeighbors = MapEmptyHood(hood);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(1/1440 * MinToHour, TIMESTEP_AGENT)) {

                int chosenCell = hood[G.rn.Int(emptyNeighbors)];

                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = naiveTcell;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                child.pd_l1 = child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
            }
            for (int run = 0; run < 3; run++) {

                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots

                // Killing logic
                for (int j = 0; j < options; j++) {
                    if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM && G.GetAgent(movdivHood[j]).bcmaExpression > 0) {
                        //G.GetAgent(movdivHood[j]).Tcell_Kill();
                        this.pd_l1 = G.boundedGaussian(10, 1, 1, 20);
                        this.type = activeTcell;
                    }
                }

                // Movement logic - only execute if no myeloma cell was encountered
                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == activeTcell) {

            if(G.TGFB.Get(Isq()) >= G.TGFBthresh){
                this.pd_1+=1;
            }

            this.myeloma_bound = false;


            int[] hood = MooreHood(true);
            // Age update every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }

            // Death from aging
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }

            // Exhaustion check
            if (this.pd_1 > pd_l1) {
                this.type = EXHT_CELL;
                return;
            }

            // -----------------------------
            // Proliferation
            // -----------------------------
            int emptyNeighbors = MapEmptyHood(hood);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) ) {

                int chosenCell = hood[G.rn.Int(emptyNeighbors)];

                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = activeTcell;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                child.pd_l1 = G.boundedGaussian(20, 1, 10, 22);
            }

            for (int run = 0; run < 3; run++) {
                // -----------------------------
                // Check neighbors for myeloma
                // -----------------------------
                int occupied = MapOccupiedHood(hood);
                simpleBoneCell target = null;

                for (int j = 0; j < occupied; j++) {

                    simpleBoneCell neighbor = G.GetAgent(hood[j]);

                    if (neighbor.type == MM && neighbor.bcmaExpression > 0) {
                        target = neighbor;
                        break;
                    }
                }
                // -----------------------------
                // Tumor interaction
                // -----------------------------
                if (target != null) {
                    this.myeloma_bound = true;
                    double kill_prob = 0;
                    if (target.mhcIExpression > 0) {
                        kill_prob = p_kill + (p_kill * 0.1);
                    }
                    if (G.rn.Double() < kill_prob) {
                        target.Tcell_Kill();
                        if (G.rn.Double() < 0.25) {
                            this.pd_1 += 1;
                        }
                        break;
                    }
                }
                // -----------------------------
                // Movement
                // -----------------------------
                int moveTo = seekCXCL9();
                if (moveTo >= 0 && G.GetAgent(moveTo) == null) {
                    MoveSQ(moveTo);
                }
            }
        }

        if (type == EXHT_CELL) {
            double EXH_RECOG_NOISE = 0.5;

            int[] hood = MooreHood(true);

            // Age update every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }

            // Age-related death
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();

                return;
            }

            if (this.pd_1 >= this.pd_l1) {
                this.Dispose();
                return;
            }

            // -----------------------------
            // Proliferation (reduced)
            // -----------------------------
            int emptyNeighbors = MapEmptyHood(hood);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(1/1440 *MinToHour, TIMESTEP_AGENT)) {

                int chosenCell = hood[G.rn.Int(emptyNeighbors)];

                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = EXHT_CELL;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(10, 1, 10, 34);
                child.pd_l1 = child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
            }
            for (int run = 0; run < 3; run++) {
                // -----------------------------
                // Check neighbors for tumor
                // -----------------------------
                int occupied = MapOccupiedHood(hood);
                simpleBoneCell target = null;

                for (int j = 0; j < occupied; j++) {

                    simpleBoneCell neighbor = G.GetAgent(hood[j]);

                    if (neighbor.type == MM && neighbor.bcmaExpression > 0) {

                        target = neighbor;
                        break;
                    }
                }

                // -----------------------------
                // Tumor interaction
                // -----------------------------
                if (target != null) {
                    // Recognition noise
                    if (G.rn.Double() >= EXH_RECOG_NOISE) {

                        double killProb = extp_kill;

                        if (G.rn.Double() < killProb) {
                            target.Tcell_Kill();
                        }
                        this.pd_1 += 1;
                    }

                }

                // -----------------------------
                // Movement (dysfunctional chemotaxis)
                // -----------------------------
                int moveTo = dysseekCXCL9();

                if (moveTo >= 0 && G.GetAgent(moveTo) == null) {
                    MoveSQ(moveTo);
                }
            }

        }

        if (type == supressorTcell) {

            int[] hood = MooreHood(true);

            // Age update every 24 hours
            if (timeStep % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }

            // Age-related death
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }

            // Self exhaustion
//            if (this.pd_1 >= this.pd_l1) {
//                this.Dispose();
//                return;
//            }
            // -----------------------------
            // Check neighbors for T cells
            // -----------------------------
            int occupied = MapOccupiedHood(hood);
            simpleBoneCell target = null;
            boolean divCheck = true;

            for (int j = 0; j < occupied; j++) {

                simpleBoneCell neighbor = G.GetAgent(hood[j]);

                if (neighbor.type == activeTcell ||
                        neighbor.type == EXHT_CELL ||
                        neighbor.type == naiveTcell) {

                    target = neighbor;
                if (neighbor.type== MM){
                    divCheck = false;
                };
                    break;
                }
            }

            // -----------------------------
            // Interaction with nearby T cells
            // -----------------------------
            if (target != null) {
                TGFB_on = true;

            }
            else {
                TGFB_on = false;
            }
            // Proliferation near T cells
            int emptyNeighbors = MapEmptyHood(hood);
            //&& Cell_Counts[12] <= 30
            if (emptyNeighbors > 0 && G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT)&& (Cell_Counts[12]) <= 50) {
                int chosenCell = hood[G.rn.Int(emptyNeighbors)];
                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = supressorTcell;
                child.tcellAge = 0;
                child.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
            }



            // -----------------------------
            // Movement toward CXCL9
            // -----------------------------
//            int options = MapEmptyHood(hood); // number of empty squares
//
//            if (options > 0) {
//                int moveTo = hood[G.rn.Int(options)]; // pick one randomly using HAL's Rand
//                MoveSQ(moveTo);
//            }
            if ((Cell_Counts[7] + Cell_Counts[8] + Cell_Counts[9]) <= 1500){
                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }

            }
            else{
                int moveToIndex = dysseekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }




        if (type == bloodVessel){
        }
    }
}
