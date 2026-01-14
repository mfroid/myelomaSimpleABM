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
    public double pmutate = 0.0;
    public double antigenLoss = Math.pow(10, -3);
    public double T_CELL_DIV_RATE = 1.0 / 1440 * (MinToHour); // T_CELL DIVISION RATE
    double CXCL9_productionRate = (2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE)/8;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE; //changed from 2.61e-10
    double CXCL9_decayRate = -.1*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE;
    double CXCL9_DiffCoef = 2700.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);//*(TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);
    double maxCXCL9 = (1.7e-9)/8;
    double IFNG_productionRate = (2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE)/8;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE; //changed from 2.61e-10
    double IFNG_decayRate = -0.005*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE;
    double IFNG_DiffCoef = 2500.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);//*(TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);
    double maxIFNG = (1.7e-9)/8;

    public static int daysPassed = 0;
    //MODEL TESTS
    public double MarrowArea;
    double convert_to_days = (MinToHour*TIMESTEP_AGENT)/(60.0*24.0); //1 ts = 6 min = 1/240 day
    int count_BA = 0;
    int init_BA = 0;
    int Nts = (int) ((numSteps)/(MinToHour*TIMESTEP_AGENT));

    public Rand rn;
    public PDEGrid2D CXCL9;
    public PDEGrid2D IFNG;

    public int[] tmoveHood = MooreHood(true);
    public ArrayList<Integer> InitBoneList = new ArrayList<>();
    public ArrayList<simpleBoneCell> AllBoneList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<simpleBoneCell> LiningList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs

    List<Object> cellSIMindex = new ArrayList<>();

    //public int eOpt;
    FileIO out;

    FileIO clones;

    FileIO bones;

    FileIO locations;

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

        //Create 2D PDE Grid for RANKL and boundary condition
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
        UIGrid IFNG_vis = new UIGrid(xDim, yDim, 2);
        UILabel days = new UILabel("days:______________________");

        if (!HEADLESS) {
            win.AddCol(0, new UILabel("Cells"));
            win.AddCol(1, days);
            win.AddCol(0, Cell_vis);
            win.AddCol(2, new UILabel("CXCL9"));
            win.AddCol(2, CXCL9_vis);
            win.AddCol(2, new UILabel("IFN- \u03B3"));
            win.AddCol(2, IFNG_vis);
            win.RunGui();
        }

        // -------------------------
        // GIF makers
        // -------------------------
        GifMaker gm_Cell_vis = new GifMaker(subfolder.concat("/CellVid.gif"), 100, true);
        GifMaker gm_CXCL9_vis = new GifMaker(subfolder.concat("/CXCL9.gif"), 100, true);
        GifMaker gm_IFNG_vis = new GifMaker(subfolder.concat("/IFNG.gif"), 100, true);

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

        double reducedMeanFraction = 0.23;
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
        double tce_reducedMeanFraction = 0.23;
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
        double no_tce_reducedMeanFraction = 0.23;
        int no_tce_totalTcells = (int) Math.floor(condition[1]*.28);
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
            else if (BIWEEKLY) {
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

            // place vessels once (unchanged)
            if (!initial_vessel) {
                int vesselNumber = 200;
                int k = 0;
                while (k < vesselNumber) {
                    int xinit = g.rn.Int(xDim);
                    int yinit = g.rn.Int(yDim);
                    if (g.GetAgent(xinit, yinit) != null && g.GetAgent(xinit, yinit).type != LINING && g.GetAgent(xinit, yinit).type != BONE) {
                        g.GetAgent(xinit, yinit).Dispose();
                    }
                    while (g.PopAt(xinit, yinit) > 0) {
                        xinit = g.rn.Int(xDim);
                        yinit = g.rn.Int(yDim);
                    }
                    simpleBoneCell c = g.NewAgentSQ(xinit, yinit);
                    vessLocations.add(c.Isq());
                    c.type = bloodVessel;
                    k++;
                }
                initial_vessel = true;
            }

            if (TCELL && initial_recruitment == false) {
                Random rand = new Random();
                int InitTcells = activeTcells; // Number of T-cells to place //160 for TCE // 60
                int InitTregs = activeTregs; // Number of suppressor Tcells to place
                int j = 0; // Counter for successfully placed T-cells
                int k = 0; // Counter for successfully placed T-cells
                int maxAttempts = 1000; // Prevent infinite loops
                int attempts = 0; // Counter for total attempts
                Collections.shuffle(vessLocations, rand);
                int[] boundaries = g.BoundaryIs();

                while (j < InitTcells && attempts < maxAttempts) {
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
                                c.type = naiveTcell;
                                if (g.rn.Double() < reducedMeanFraction) {
                                    c.type = EXHT_CELL;
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                } else {
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                }
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
                                c.type = naiveTcell;
                                if (g.rn.Double() < reducedMeanFraction) {
                                    c.type = EXHT_CELL;
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                } else {
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                }
                                j++;
                                break;
                            }
                            tries++;
                        }
                    }

                    if (attempts >= maxAttempts) {
                    }
                }

                while (k < InitTregs && attempts < maxAttempts) {
                    attempts++;

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
                                c.type = supressorTcell;
                                c.pd_1 = g.boundedGaussian(10, 1, 10, 20);
                                k++;
                                placed = true;
                            }
                        }
                    }

                    // Fallback: pick random boundary if no vessel-adjacent space was found
                    if (!placed) {
                        int tries = 0;
                        while (tries < boundaries.length) {
                            int randomBoundaryIdx = rand.nextInt(boundaries.length);
                            int candidate = boundaries[randomBoundaryIdx];
                            if (g.GetAgent(candidate) == null) {
                                simpleBoneCell c = g.NewAgentSQ(candidate);
                                c.type = supressorTcell;
                                c.pd_1 = g.boundedGaussian(10, 1, 10, 20);
                                k++;
                                break;
                            }
                            tries++;
                        }
                    }

                }

                // Update flags to indicate T-cell recruitment is complete
                initial_recruitment = true;
            }

            if (TCELL && initial_recruitment && daysPassed > 0) {
                int InitTcells;
                if (Cell_Counts[11]+Cell_Counts[12] < 100) {
                    InitTcells = (int) Math.floor(g.boundedGaussian(10, 1, 0, 10));
                } else {
                    InitTcells = (int) Math.floor(g.boundedGaussian(1, 1, 0, 10));

                }



                int j = 0;
                int maxAttempts = 1000;
                int attempts = 0;
                Random rand = new Random();
                Collections.shuffle(vessLocations, rand);
                int num_placed = 0;
                int[] boundaries = g.BoundaryIs();


                while (j < InitTcells && attempts < maxAttempts) {
                    attempts++;

                    boolean placed = false;

                    // Select a random blood vessel location
                    int randomVessLocation = vessLocations.get(rand.nextInt(vessLocations.size()));

                    if (g.GetAgent(randomVessLocation).type == bloodVessel) {
                        int[] movdivHood = MooreHood(true);
                        int emptyNeighbors = g.MapEmptyHood(movdivHood, randomVessLocation);

                        if (emptyNeighbors > 0) {
                            int chosenIndex = g.rn.Int(emptyNeighbors);
                            int chosenCell = movdivHood[chosenIndex];

                            if (g.GetAgent(chosenCell) == null) {
                                simpleBoneCell c = g.NewAgentSQ(chosenCell);
                                if (g.rn.Double() <.2){
                                    c.type = EXHT_CELL;
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                    j++;
                                    placed = true;
                                    num_placed+= 1;

                                }
                                else{
                                c.type = naiveTcell;
                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                j++;
                                placed = true;
                                num_placed+= 1;
                                }
                            }
                        }
                    }

                    // Fallback to boundary if not placed
                    if (!placed) {
                        int tries = 0;
                        while (tries < boundaries.length) {
                            int randomBoundaryIdx = rand.nextInt(boundaries.length);
                            int chosenCell = boundaries[randomBoundaryIdx];
                            if (g.GetAgent(chosenCell) == null) {
                                simpleBoneCell c = g.NewAgentSQ(chosenCell);
                                if (g.rn.Double() <.2){
                                    c.type = EXHT_CELL;
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                    j++;
                                    placed = true;
                                    num_placed+= 1;

                                }
                                else{
                                    c.type = naiveTcell;
                                    c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                    c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                                    j++;
                                    placed = true;
                                    num_placed+= 1;
                                }
                            }
                            tries++;
                        }
                    }
                }

            }

            if (TCELL && initial_recruitment && daysPassed > 0) {
                int InitTcells;
                if (Cell_Counts[12] < activeTregs) {
                    InitTcells = (int) Math.floor(g.boundedGaussian((activeTregs - Cell_Counts[12]), 1, activeTregs - Cell_Counts[12]-1, activeTregs - Cell_Counts[12] + 1));
                } else {
                    //InitTcells = (int) g.boundedGaussian(1, 1, 1, 10);
                    InitTcells = 0; // as per your original logic
                }

                int j = 0;
                int maxAttempts = 1000;
                int attempts = 0;
                Random rand = new Random();
                Collections.shuffle(vessLocations, rand);
                int num_placed = 0;
                int[] boundaries = g.BoundaryIs();

                while (j < InitTcells && attempts < maxAttempts) {
                    attempts++;
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
                                c.type = supressorTcell;
                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                num_placed+=1;
                                j++;
                                placed = true;

                            }
                        }
                    }

                    // Fallback: try to place on random boundary cell
                    if (!placed) {
                        int tries = 0;
                        while (tries < boundaries.length) {
                            int randomBoundaryIdx = rand.nextInt(boundaries.length);
                            int chosenCell = boundaries[randomBoundaryIdx];
                            if (g.GetAgent(chosenCell) == null) {
                                simpleBoneCell c = g.NewAgentSQ(chosenCell);
                                c.type = supressorTcell;
                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                num_placed+=1;
                                j++;
                                break;
                            }
                            tries++;
                        }
                    }

                }

            }

            if (g.rn.Double() < dailyTregIncreasenProb) {
                if (activeTregs < maxTregs) {
                    activeTregs += 1;
                }
            }
            // Model step and drawing (same as before)
            g.ModelStep(i, Cell_Counts, simID);

            if (!HEADLESS) {
                g.Draw(Cell_vis, days, i, simID);
                g.DrawCXCL9(CXCL9_vis);
                g.DrawIFNG(IFNG_vis);
            }

            // daily recording: your original condition:
            if (i % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                g.RecordOut(g.out, i, TREATMENT_ON, BORTEZOMIB, MYELOMA);
                g.RecordClones(g.clones, i);
                g.RecordBones(g.bones, i);
                g.RecordLocs(g.locations, i);
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
        gm_IFNG_vis.Close();

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

        paramsOut = new FileIO(projPath + "paramsOut.txt",mode);


        if(mode=="w") {

            out.Write("Timestep" + "," + "BONE" + "," + "pOB" + "," + "aOB" + "," + "pOC" + "," + "aOC" + "," + "MSC" + "," + "LINING" + "," + "S_MM" + "," + "R_MM"+"," +"AL_MM"+ ","+"TCell"+","+"ExtTcell"+"," +"T-reg"+"," + "Naive Tcell" +","+ "TREATMENT_ON" + "," + "BORTEZOMIB" + "," + "MYELOMA" + "\n");
            clones.Write("Timestep" + "," + "SimID" + "," + "MHCI" + "," + "BCMA" + "\n");
            bones.Write("Timestep" + "," + "SimID" + "," + "Position" + "\n");
            locations.Write("Timestep" + "," + "SimID" + "," + "Position" + ","+"Type"+"\n");
            paramsOut.Write("Header");
        }

    }

    public void closeFileIO () {

        out.Close();
        clones.Close();
        bones.Close();
        locations.Close();
        paramsOut.Close();
//        params.Close();
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


        Queue<int[]> cellQueue = new LinkedList<>(); // Queue to manage cluster growth
        Set<String> visited = new HashSet<>(); // Track visited cells to prevent duplicates

        // Step 1: Seed the initial myeloma cell near the bone
        boolean seedPlaced = false;
        while (!seedPlaced) {
            int xInit = rn.Int(xDim);
            int yInit = rn.Int(yDim);

            // Check if the location is near a bone cell within the specified proximity
            boolean isNearBone = false;
            for (int xi = Math.max(0, xInit - boneProximityDistance); xi <= Math.min(xDim - 1, xInit + boneProximityDistance); xi++) {
                for (int yi = Math.max(0, yInit - boneProximityDistance); yi <= Math.min(yDim - 1, yInit + boneProximityDistance); yi++) {
                    if (GetAgent(xi, yi) != null && GetAgent(xi, yi).type == BONE) {
                        double distance = Math.sqrt(Math.pow(xi - xInit, 2) + Math.pow(yi - yInit, 2));
                        if (distance <= boneProximityDistance) {
                            isNearBone = true;
                            break;
                        }
                    }
                }
                if (isNearBone) break;
            }

            // Place the initial myeloma cell
            if (isNearBone && PopAt(xInit, yInit) == 0) {
                NewAgentSQ(xInit, yInit).type = MM; // Seed the initial myeloma cell
                GetAgent(xInit, yInit).bcmaExpression = 1;
                cellQueue.add(new int[]{xInit, yInit}); // Add to queue for cluster growth
                visited.add(xInit + "," + yInit); // Mark as visited
                placedMyelomaCells++;
                seedPlaced = true;
            }
        }

        // Step 2: Grow the cluster using a true circular expansion
        while (placedMyelomaCells < myelomaCellsToPlace && !cellQueue.isEmpty()) {
            int[] currentCell = cellQueue.poll(); // Get the next cell from the queue
            int xCurrent = currentCell[0];
            int yCurrent = currentCell[1];

            // Randomly sample points within a circular radius
            for (int i = 0; i < 8; i++) { // Limit to 8 random points per cell to keep placement organic
                double angle = rn.Double() * 2 * Math.PI; // Random angle
                double radius = rn.Double() * 2.0; // Random radius (adjust scale for tighter/looser packing)
                int xNeighbor = xCurrent + (int) Math.round(radius * Math.cos(angle));
                int yNeighbor = yCurrent + (int) Math.round(radius * Math.sin(angle));

                // Check if the neighbor is within bounds and unvisited
                if (xNeighbor >= 0 && xNeighbor < xDim && yNeighbor >= 0 && yNeighbor < yDim &&
                        !visited.contains(xNeighbor + "," + yNeighbor)) {

                    // Place the myeloma cell if the location is unoccupied
                    if (PopAt(xNeighbor, yNeighbor) == 0) {
                        NewAgentSQ(xNeighbor, yNeighbor).type = MM; // Place the cell
                        GetAgent(xNeighbor, yNeighbor).bcmaExpression = 1;
                        if (rn.Double() < bcmaNegfraction){
                            GetAgent(xNeighbor, yNeighbor).bcmaLoss = true;
                            GetAgent(xNeighbor, yNeighbor).bcmaExpression = rn.Double();


                        }
                        cellQueue.add(new int[]{xNeighbor, yNeighbor}); // Add to the queue
                        visited.add(xNeighbor + "," + yNeighbor); // Mark as visited
                        placedMyelomaCells++;

                        // Stop if we've placed all required cells
                        if (placedMyelomaCells >= myelomaCellsToPlace) {
                            break;
                        }
                    }
                }
            }
        }



    }

    public void ModelStep(int time, double [] Cell_Counts, int simID) {

        //STEP 0: UPDATE GRIDTICK
        /////////////////////////////////////////////////
        //STEP 1: REACTION-DIFFUSION EQUATION FOR RANKL//
        /////////////////////////////////////////////////
        int i=0;
        double stol = 1.0e-6;//1e-6; //steady-state tolerance


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

        for (int x = 0; x < IFNG.xDim; x++) {
            for (int y = 0; y < IFNG.yDim; y++) {
                if (GetAgent(x,y)!=null) {
                    if (GetAgent(x, y).myeloma_bound==true || GetAgent(x, y).type == naiveTcell){
                        IFNG.Add(x, y, IFNG_productionRate / maxIFNG);
                    }
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    IFNG.Set(x, y, 0.1 * IFNG_DiffCoef);
                } else if (x!=xDim-1 && GetAgent(x+1, y) != null && (GetAgent(x+1, y).type == BONE || GetAgent(x+1, y).type == LINING)){
                    IFNG.Set(x, y, 0.1 * IFNG_DiffCoef);
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    IFNG.Set(x, y, 0.1 * IFNG_DiffCoef);
                } else if (y!=yDim-1 && GetAgent(x, y+1) != null && (GetAgent(x, y+1).type == BONE || GetAgent(x, y+1).type == LINING)){
                    IFNG.Set(x, y, 0.1 * IFNG_DiffCoef);
                }
            }
        }

        //IFNG Diffusion
        IFNG.DiffusionADI(IFNG_DiffCoef);

        //Natural Decay of IFNG
        IFNG.MulAll(IFNG_decayRate);
        IFNG.Update();


        //System.out.println("max RANKL "+rmax);
        //System.out.println("max TGFB "+tmax);


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
                double mhcI = c.mhcIExpression;
                double bcma = c.bcmaExpression;
                writeHere.Write(time + "," + simID + "," + mhcI + ","+ bcma + "\n");
            }
        }
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
        double[] cts = new double[13];

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
                {3000, 208}
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

    public int seekIFNG() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] IFNG_levels = new double[9]; // stores IFNG levels at the 8 directions and the center

        // Extracting IFNG levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            IFNG_levels[i] = G.IFNG.Get(G.tmoveHood[i]);
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
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[1] - IFNG_levels[2]);
                        break;
                    case 2: // left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[1] - IFNG_levels[2]);
                        break;
                    case 3: // up
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[3] - IFNG_levels[4]);
                        break;
                    case 4: // down
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[3] - IFNG_levels[4]);
                        break;
                    case 5: // top-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[5] - IFNG_levels[6]);
                        break;
                    case 6: // top-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[5] - IFNG_levels[6]);
                        break;
                    case 7: // bottom-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[7] - IFNG_levels[8]);
                        break;
                    case 8: // bottom-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxIFNG) / 8 * (IFNG_levels[7] - IFNG_levels[8]);
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

    void Init() {
        if (type == BONE || type==LINING) { //Include aOB because they will turn into bone
            G.count_BA++;
        }
    }


    public boolean MarrowInHood () {
        int[] MarrowHood = VonNeumannHood(false); //4 neighbors
        int options = MapEmptyHood(MarrowHood); //occupied positions surrounding pOB
        return options > 0;
    }

    double tcellAge=0;
    double lifespan=0;
    double currentAge=0;
    double pd_1 = 0;
    double pd_l1 = 0;
    void Tcell_Kill() {
        if (type == MM) {
            this.Dispose();}
    }
    boolean myeloma_bound = false;
    boolean mhcLoss = false;
    double bcmaExpression = 1.0;
    double mhcIExpression = 1.0;
    int simulationID;

    public void CellStep(int time, double [] Cell_Counts, int simID) {


        ///////////
        //MYELOMA//
        ///////////
        if (type == MM) {

            if (time ==0){

            this.simulationID = simID;
            }
//            if (this.simulationID!=simID) {
//                color = RGB(255, 255, 0);
//            }
            double rn_BirthDeath = G.rn.Double();
            double pdiv;
            double pdeath = G.MM_DEATH;

            double scaleFactor = 0.7 + (0.3 * this.bcmaExpression); // bcma=0 → 0.7, bcma=1 → 1.0
            pdiv = ProbScale((1.0 / 2400) * scaleFactor, MinToHour);
            ///////////
            //MM dies//
            ///////////
            if (rn_BirthDeath < ProbScale(pdeath, TIMESTEP_AGENT)) {//rn_BirthDeath < pdeath){
                //MM dies
                color = WHITE;
                Dispose();
            }
            else if (rn_BirthDeath < ProbScale(pdiv, TIMESTEP_AGENT)) {

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

                        if (G.rn.Double() < G.antigenLoss){
                            //double[] possibleValues = {0.0, 0.25, 0.5, 0.75};
                            //int index = G.rn.Int(possibleValues.length); // random integer 0–3
                            //double newExpression = possibleValues[index];
                            child.mhcIExpression = 0;
                            //child.mhcLoss = true;
                        }
                        else{child.mhcIExpression = this.mhcIExpression;
                        }

                        // Antigen loss mutation
                        if (G.rn.Double() < G.antigenLoss && bcmaLoss==false) {
                            // Pick one of four discrete expression levels: 0, 0.25, 0.5, or 0.75
                            //double[] possibleValues = {0.0, 0.25, 0.5, 0.75};
                            //int index = G.rn.Int(possibleValues.length); // random integer 0–3
                            //double newExpression = possibleValues[index];
                            if (child != null) {
                                child.bcmaExpression = 0;
                                //child.bcmaLoss = true;
                            }
                        } else {
                            child.bcmaExpression = this.bcmaExpression;
                            child.mhcLoss = mhcLoss;
                        }
                    }
                }
            }

            int[] Hood = MooreHood(true); // For division and movement
            int options = MapOccupiedHood(Hood); // Mapping occupied spots
            double x = 4.73e-8;
            for (int j = 0; j < options; j++) {
                if (G.GetAgent(Hood[j]) != null && G.GetAgent(Hood[j]).type == BONE || G.GetAgent(Hood[j]).type == LINING) {
                    if(G.rn.Double() < x){
                    G.GetAgent(Hood[j]).Dispose();
                        break;
                    }

                }
            }
        }

        if (type == naiveTcell) {
            this.myeloma_bound = false;
            boolean encounteredMyeloma = false; // Track if a myeloma cell is encountered

            // Increment T cell age every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }

            for (int run = 0; run < 6; run++) {
                // Transition to exhausted T cell if pd_l1 exceeds pd_1
                if (this.pd_l1 >= this.pd_1) {
                    this.type = EXHT_CELL;
                    break;
                }
                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots

                // Age-related death logic
                if (this.tcellAge >= this.lifespan) {
                    this.Dispose();
                    break;
                }

                else {
                    double x = ProbScale(.08, TIMESTEP_AGENT);
                    //double x = .12;4

                    // Killing logic
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM && G.GetAgent(movdivHood[j]).bcmaExpression>0) {
                            encounteredMyeloma = true; // Stop movement if any myeloma cell is found
                            this.myeloma_bound = true;
                            if (G.GetAgent(movdivHood[j]).type == MM && (G.GetAgent(movdivHood[j]).mhcIExpression>0)){
                                this.type = activeTcell;
                                break;
                            }
                            else if (G.GetAgent(movdivHood[j]).bcmaExpression>0 && (G.GetAgent(movdivHood[j]).mhcIExpression==0) ) {
                                if (G.rn.Double() < x) {
                                    G.GetAgent(movdivHood[j]).Tcell_Kill();
                                    this.pd_l1 += 1;
                                    this.type = EXHT_CELL;
                                }
                                break; // Exit the inner loop as soon as a myeloma cell is encountered
                            }
                        }
                    }
                }

                // If a myeloma cell is encountered, break out of the outer loop immediately
                if (encounteredMyeloma) {
                    break;
                }

                // Movement logic - only execute if no myeloma cell was encountered
                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == activeTcell) {
            this.myeloma_bound = false;
            boolean encounteredMyeloma = false; // Track if a myeloma cell is encountered

            // Increment T cell age every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }

            for (int run = 0; run < 6; run++) {
                // Transition to exhausted T cell if pd_l1 exceeds pd_1
                if (this.pd_l1 >= this.pd_1) {
                    this.type = EXHT_CELL;
                    break;
                }
                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots

                // Age-related death logic
                if (this.tcellAge >= this.lifespan) {
                    this.Dispose();
                    break;
                }

                else {
                    double x = ProbScale(.12, TIMESTEP_AGENT);
                    //double x = .12;4

                    // Killing logic
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM && (G.GetAgent(movdivHood[j]).mhcIExpression>0 || G.GetAgent(movdivHood[j]).bcmaExpression>0)) {
                            encounteredMyeloma = true; // Stop movement if any myeloma cell is found
                            this.myeloma_bound = true;
                            if (G.rn.Double() < x) {
                                G.GetAgent(movdivHood[j]).Tcell_Kill();
                                this.pd_l1 += 1;
                            }
                            int emptyNeighbors = MapEmptyHood(MooreHood(true)); // mapping empty spots
                            if (emptyNeighbors > 0 &&  G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) ) {
                                for (int i = 0; i < emptyNeighbors; i++) {
                                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                                    int chosenCell = MooreHood(true)[chosenIndex]; // Get the chosen empty cell
                                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                        simpleBoneCell child = G.NewAgentSQ(chosenCell);
                                        child.type = this.type;
                                        child.tcellAge = 0;
                                        child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                                        child.pd_l1 = child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
                                        break;
                                    }
                                }
                            }
                            break; // Exit the inner loop as soon as a myeloma cell is encountered
                        }
                    }
                }

                // If a myeloma cell is encountered, break out of the outer loop immediately
                if (encounteredMyeloma) {
                    break;
                }

                // Movement logic - only execute if no myeloma cell was encountered
                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }


        if (type == EXHT_CELL) {
            this.myeloma_bound = false;
            boolean encounteredMyeloma = false; // Track if a myeloma cell is encountered

            // Increment T cell age every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }
//            int emptyNeighbors = MapEmptyHood(MooreHood(true)); // mapping empty spots
//            if (emptyNeighbors > 0 &&  G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) ) {
//                for (int i = 0; i < emptyNeighbors; i++) {
//                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
//                    int chosenCell = MooreHood(true)[chosenIndex]; // Get the chosen empty cell
//                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
//                        simpleBoneCell child = G.NewAgentSQ(chosenCell);
//                        child.type = this.type;
//                        child.tcellAge = 0;
//                        child.pd_l1 = child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
//                        break;
//                    }
//                }
//            }


            for (int run = 0; run < 6; run++) {
                // Transition to exhausted T cell if pd_l1 exceeds pd_1
                if (this.pd_l1 >= this.pd_1) {
                    this.Dispose();
                    break;
                }
                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots

                // Age-related death logic
                if (this.tcellAge >= this.lifespan) {
                    this.Dispose();
                    break;
                }

                else {

                    double x = ProbScale( 0.12*(.5),TIMESTEP_AGENT);
                    //double x = .12/2;
                    // Killing logic
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM && (G.GetAgent(movdivHood[j]).mhcIExpression>0 || G.GetAgent(movdivHood[j]).bcmaExpression>0)) {
                            encounteredMyeloma = true; // Stop movement if any myeloma cell is found
                            this.myeloma_bound = true;
                            if (G.rn.Double() < x) {
                                this.pd_l1 += 1;
                                G.GetAgent(movdivHood[j]).Tcell_Kill();
                            }
                            break; // Exit the inner loop as soon as a myeloma cell is encountered
                        }
                    }
                }

                // If a myeloma cell is encountered, break out of the outer loop immediately
                if (encounteredMyeloma) {
                    break;
                }

                // Movement logic - only execute if no myeloma cell was encountered
                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }


        if (type == supressorTcell) {
            this.lifespan = G.boundedGaussian(40, 1, 40, 44);

            if (timeStep % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }



            for (int run = 0; run < 6; run++) {

                int[] movdivHood = MooreHood(true); // for division and movement
                int options = MapOccupiedHood(movdivHood); // mapping occupied spots


                // Age-related death logic
                if (this.tcellAge >= this.lifespan) {
                    this.Dispose();
                    break;
                }
                if (this.pd_l1 >= this.pd_1) {
                    this.Dispose();
                    break;
                }
                else {
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && (G.GetAgent(movdivHood[j]).type == activeTcell || G.GetAgent(movdivHood[j]).type == EXHT_CELL|| G.GetAgent(movdivHood[j]).type == naiveTcell) && G.rn.Double() < G.boundedGaussian(.9, .1, 0, 1)) {
                            int emptyNeighbors = MapEmptyHood(MooreHood(true)); // mapping empty spots
                            if (emptyNeighbors > 0 &&  G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT)) {
                                for (int i = 0; i < emptyNeighbors; i++) {
                                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                                    int chosenCell = MooreHood(true)[chosenIndex]; // Get the chosen empty cell
                                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                        simpleBoneCell child = G.NewAgentSQ(chosenCell);
                                        child.type = this.type;
                                        child.tcellAge = 0;
                                        child.pd_1 = G.boundedGaussian(10, 1, 10, 20);
                                        break;
                                    }
                                }
                            }
                            G.GetAgent(movdivHood[j]).pd_l1 = G.GetAgent(movdivHood[j]).pd_l1+1;
                            this. pd_l1 +=1;
                            break; // Break out of for loop
                        }
                    }

                    // Movement logic only executes if no T cell was killed
                    int moveToIndex = seekIFNG();
                    if (G.GetAgent(moveToIndex) == null) {
                        MoveSQ(moveToIndex);
                    }
                }
            }
        }

        if (type == bloodVessel){
        }
    }
}
