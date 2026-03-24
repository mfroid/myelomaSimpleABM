package simpleMyelomaABM;
import HAL.GridsAndAgents.*;
import HAL.Gui.*;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;;
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
    public static boolean TCELL = true;
    public static boolean TCE = true;
    public static boolean BIWEEKLY = false;
    public static boolean EMDR = false;

    public static boolean TREATMENT_ON = false; //this is to control treatment on/off timer in MAIN

    //CLUSTER
    public static boolean PARAM_SWEEP = false; //use when importing parameters to loop through
    public static boolean runPar = true;
    public static boolean HEADLESS = false; //use true with cluster
    public static boolean LOCAL = true; // use false with cluster
    public static double numSteps = 1.0*365.0*24.0*60.0; // years the model will run
    public static int numSims = 10; //Number of Simulations
    public final static int BONE = RGB256(255,255,250),
            LINING = RGB256(64,106,151),
            MM = RGB256(0,128,0),
            activeTcell = RGB256(17, 150, 150),
            EXHT_CELL=RGB256(200, 50, 250),
            supressorTcell =RGB256(255, 165, 0),
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
    double CXCL9_decayRate = -.2*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//*(TIMESTEP_AGENT)/N_TIMESTEP_PDE;
    double CXCL9_DiffCoef = 2700.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);//*(TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE);
    double maxCXCL9 = (2.07e-9)/8;

    double PERF_productionRate = 2.04e-9 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE; // 2.04e-9 //dts; PERF production rate
    static double PERF_basalRate = 2.04e-11 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE; // 2.04e-11 //dts; basal PERF production rate (not by OC)
    static double PERF_decayRate = -0.35 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE; // dts; BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    double PERF_DiffCoef = 780.0 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP * N_TIMESTEP_PDE); // dts
    double Extra_PERFtime = 4320.0 / (MinToHour * TIMESTEP_AGENT); // 4320 min = 720 ts = 3 days
    static double maxPERF = 8.7e-10; // 1.12e-9; //1.78e-9;//1.21e-9; //6.4e-10

    //Cytokines
    double TGFB_productionRate = 2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//2.04e-9//dts; RANKL production rate
    static double TGFB_basalRate = 2.04e-11*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //2.04e-11//dts; basal TGFB production rate (not by OC)
    static double TGFB_decayRate = -0.35*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //dts; BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    double TGFB_DiffCoef = 780.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE); //dts
    double Extra_TGFBtime = 4320.0/(MinToHour*TIMESTEP_AGENT); //4320 min = 720 ts = 3 days
    static double maxTGFB = 8.7e-10; //1.12e-9; //1.78e-9;//1.21e-9; //6.4e-10

    double Ts = TGFB_basalRate/(Math.abs(TGFB_decayRate)*maxTGFB); //Basal TGFB
    double TGFBthresh = (1.05)* TGFB_basalRate/(Math.abs(TGFB_decayRate)* maxTGFB);//was 0.01; now 5% increase from basal
    double MaxTdiff; //max relative change in TGFB
    double tmax;//=0;
    boolean [][] extraTGFB = new boolean[xDim][yDim];
    double [][] TGFBtimer = new double[xDim][yDim];

    //INHIBITORS AND TREATMENT
    double dose = 1.0; //1.0 ;//1.0; //1.0
    public int Tx_Interval = (int) (30240.0/(MinToHour*TIMESTEP_AGENT)); //(5760.0/(TIMESTEP_AGENT)); //5760 = 4 day; //30240 min = 21 days
    public int Tx_Duration = (int) (30240.0/(MinToHour*TIMESTEP_AGENT));//(20160.0/(TIMESTEP_AGENT)); //1440 = 1 day; 4320 = 3 days; 525600 = 365 days; 20160 min = 14 days
    public static int daysPassed = 0;
    //MODEL TESTS
    public double MarrowArea;
    double convert_to_days = (MinToHour*TIMESTEP_AGENT)/(60.0*24.0); //1 ts = 6 min = 1/240 day
    int count_BA = 0;
    int init_BA = 0;
    int Nts = (int) ((numSteps)/(MinToHour*TIMESTEP_AGENT));

    public Rand rn;
    public PDEGrid2D CXCL9;
    public PDEGrid2D TGFB;

    public PDEGrid2D PERF;

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
        TGFB = new PDEGrid2D(xDim, yDim,true,true);
        PERF = new PDEGrid2D(xDim, yDim,true,true);
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

    // --- Add snapshot export helper in simpleBoneGrid (instance method) ---
    public CellSnapshot exportSnapshot(simpleBoneCell c) {
        // Capture the fields you care about. Extend as needed.
        return new CellSnapshot(
                c.type,
                c.bcmaLoss,
                c.simulationID
        );
    }

    // --- Add snapshot import helper in simpleBoneGrid (instance method) ---
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
        UIGrid TGFB_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
        UIGrid PERF_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
        UILabel days = new UILabel("days:______________________");

        if (!HEADLESS) {
            win.AddCol(0, new UILabel("Cells"));
            win.AddCol(1, days);
            win.AddCol(0, Cell_vis);
            win.AddCol(2, new UILabel("CXCL9"));
            win.AddCol(2, CXCL9_vis);
            win.AddCol(2, new UILabel("TGF-\u03B2"));
            win.AddCol(2, TGFB_vis);
            win.AddCol(2, new UILabel("Perforin"));
            win.AddCol(2, PERF_vis);
            win.RunGui();
        }

        // -------------------------
        // GIF makers
        // -------------------------
        GifMaker gm_Cell_vis = new GifMaker(subfolder.concat("/CellVid.gif"), 100, true);
        GifMaker gm_CXCL9_vis = new GifMaker(subfolder.concat("/CXCL9.gif"), 100, true);
        GifMaker gm_TGFB_vis = new GifMaker(subfolder.concat("/TGFB.gif"), 100, true);
        GifMaker gm_PERF_vis = new GifMaker(subfolder.concat("/PERF.gif"), 100, true);
        // -------------------------
        // Bone file selection and grid creation
        // -------------------------
        String Bone_Filename = null;
        if (LOCAL) {
            Bone_Filename = "/Users/80024703/Desktop/code/Bone/BAout_2020May5_Sim14.csv";
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

        // -------------------------
        // main time loop (your original for loop)

        for (int i = 0; i < g.Nts; i++) {

            double[] Cell_Counts = g.CellCounts();
            timeStep = i;


            if (!HEADLESS) {
                win.TickPause(10);
            }

            if (TCELL && i ==0) {

                int InitTcells = 250;

//                if ( (Cell_Counts[10]+Cell_Counts[11]+Cell_Counts[12]+Cell_Counts[13]) < InitTcells){
//                    InitTcells = 250;
//                }
//                else{
//                    InitTcells = 10;
//                }
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

//                        if (g.rn.Double() < exhaustedFraction) {
//                            c.type = activeTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(1, 1, 1, 2);
//                        } else if (g.rn.Double() < .10 && i ==0){
//                            c.type = supressorTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
//                        } else{
//                            c.type = naiveTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
//                        }
                        if(g.rn.Double() < .1){
                            c.type = supressorTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        } else if (g.rn.Double() < .25){
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }else{
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
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

//                        if (g.rn.Double() < exhaustedFraction) {
//                            c.type = activeTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(1, 1, 1, 2);
//                        } else if (g.rn.Double() < .10 && i ==0) {
//                            c.type = supressorTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
//                        } else{
//                            c.type = naiveTcell;
//                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
//                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
//                        }
                        if(g.rn.Double() < .1){
                            c.type = supressorTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        } else if (g.rn.Double() < .25){
                            c.type = EXHT_CELL;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);
                        }else{
                            c.type = naiveTcell;
                            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
                            c.lifespan = g.boundedGaussian(30, 1, 30, 34);

                        }
                        j++;
                    }
                }
            }

            // compute days/weeks
            int stepsPerDay = (int)(24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT));
            int daysPassedCalc = i / stepsPerDay;
            int weeksPassed = daysPassedCalc / 7;

            // give drug once every OTHER week (week 0, 2, 4, ...)
            if (BIWEEKLY &&
                    daysPassedCalc % 7 == 0 &&      // first day of the week
                    i % stepsPerDay == 0 &&         // start of the day
                    weeksPassed % 2 == 0) {         // every other week

                TCE = true;

            } else if (BIWEEKLY) {
                TCE = false;
            }
            // Model step and drawing (same as before)
            g.ModelStep(i, Cell_Counts, simID);

            if (!HEADLESS) {
                g.Draw(Cell_vis, days, i, simID);
                g.DrawCXCL9(CXCL9_vis);
                g.DrawTGFB(TGFB_vis);
                g.DrawPERF(PERF_vis);
            }
            // daily recording: your original condition:
            if (i % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
//                System.out.println("_______________________________");
//                System.out.println(i/240);
//                System.out.println(Cell_Counts[10]+Cell_Counts[13]);
//                System.out.println("_______________________________");
                g.RecordOut(g.out, i);
                g.RecordClones(g.clones, i);
                g.RecordBones(g.bones, i);
                //g.RecordLocs(g.locations, i);

//                // Export myeloma cells for transfer
//                for (simpleBoneCell c : g) {
//                    if (c.type == MM & runPareShare) {
//                        if (g.rn.Double() < 0.01) { // tune probability
//                            CellSnapshot snap = g.exportSnapshot(c);
//                            myelomaTransferQueue.add(snap);
//                            c.Dispose(); // remove from original grid
//                        }
//                    }
//                }
            }
        } // end time loop


        // -------------------------
        // cleanup: close files, gifs, window
        // -------------------------
        g.closeFileIO();
        gm_Cell_vis.Close();
        gm_CXCL9_vis.Close();
        gm_TGFB_vis.Close();
        gm_PERF_vis.Close();

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


        if(mode=="w") {

            out.Write("Timestep" + "," + "BONE" + "," + "LINING" + "," + "S_MM" + "," + "R_MM"+"," +"AL_MM"+ ","+"TCell"+","+"ExtTcell"+"," +"T-reg" +","+ "Naive Tcell" +"\n");
            clones.Write("Timestep" + "," + "SimID" + "," + "MHCI" + "," + "BCMA" + "\n");
            bones.Write("Timestep" + "," + "SimID" + "," + "Position" + "\n");
            locations.Write("Timestep" + "," + "SimID" + "," + "Position" + ","+"Type"+"\n");
        }

    }

    public void closeFileIO () {

        out.Close();
        clones.Close();
        bones.Close();
        locations.Close();
//        params.Close();
    }

    public void SetParams(int prow, ArrayList<String> param_list){
        //returns an array list of all lines from the file as stringsftype == lining

        String[] split_param_list = param_list.get(prow).split(",");


        pmutate = Double.parseDouble(split_param_list[2]);
        EMDR = Boolean.parseBoolean(split_param_list[3]);
        dose = Double.parseDouble(split_param_list[4]);
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
        String[] split_input_data = input_data.get(0).split(",");

        //Place bone
        for (int index = 1; index < split_input_data.length; index++) {
            NewAgentSQ(Integer.parseInt(split_input_data[index])).type = BONE;
            GetAgent(Integer.parseInt(split_input_data[index])).Init();
            InitBoneList.add(Integer.parseInt(split_input_data[index]));
            AllBoneList.add(GetAgent(Integer.parseInt(split_input_data[index])));
        }
        for (int index = 1; index < split_input_data.length; index++) {
            if (GetAgent(Integer.parseInt(split_input_data[index])).MarrowInHood() == true) {
                GetAgent(Integer.parseInt(split_input_data[index])).type = LINING;
                GetAgent(Integer.parseInt(split_input_data[index])).liningAge = TURNOVER_TIME;
                LiningList.add(GetAgent(Integer.parseInt(split_input_data[index])));
            }
        }
        init_BA = InitBoneList.size();
        MarrowArea = (xDim * yDim) - init_BA;//(xDimBone*yDimBone); //0.12 Bone, 0.88 Marrow
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


            if (rn.Double() < bcmaNegfraction) {
                GetAgent(x, y).bcmaLoss = true;
                GetAgent(x, y).bcmaExpression = rn.Double();
            }

            placedMyelomaCells++;
        }




    }

    public void ModelStep(int time, double [] Cell_Counts, int simID) {

        //STEP 0: UPDATE GRIDTICK
        /////////////////////////////////////////////////
        //STEP 1: REACTION-DIFFUSION EQUATION FOR RANKL//
        /////////////////////////////////////////////////

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


        for (int x = 0; x < TGFB.xDim; x++) {
            for (int y = 0; y < TGFB.yDim; y++) {

                if (GetAgent(x,y) != null &&
                        GetAgent(x, y).type == MM &&
                        GetAgent(x, y).TGFB_on == true)  {

                    TGFB.Add(x, y, TGFB_productionRate / maxTGFB);

                }

                if (GetAgent(x, y) != null &&
                        (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {

                    TGFB.Set(x, y, 0.1 * TGFB_DiffCoef);

                } else if (x != xDim - 1 &&
                        GetAgent(x + 1, y) != null &&
                        (GetAgent(x + 1, y).type == BONE || GetAgent(x + 1, y).type == LINING)) {

                    TGFB.Set(x, y, 0.1 * TGFB_DiffCoef);
                }

                if (GetAgent(x, y) != null &&
                        (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {

                    TGFB.Set(x, y, 0.1 * TGFB_DiffCoef);

                } else if (y != yDim - 1 &&
                        GetAgent(x, y + 1) != null &&
                        (GetAgent(x, y + 1).type == BONE || GetAgent(x, y + 1).type == LINING)) {

                    TGFB.Set(x, y, 0.1 * TGFB_DiffCoef);
                }
            }
        }

        // TGFB Diffusion
        TGFB.DiffusionADI(TGFB_DiffCoef);

        // Natural Decay of TGFB
        TGFB.MulAll(TGFB_decayRate);
        TGFB.Update();
        for (int x = 0; x < PERF.xDim; x++) {
            for (int y = 0; y < PERF.yDim; y++) {

                if (GetAgent(x,y) != null &&
                        GetAgent(x, y).type == MM &&
                        GetAgent(x, y).PERF_on == true)  {

                    PERF.Add(x, y, PERF_productionRate / maxPERF);

                }

                if (GetAgent(x, y) != null &&
                        (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {

                    PERF.Set(x, y, 0.1 * PERF_DiffCoef);

                } else if (x != xDim - 1 &&
                        GetAgent(x + 1, y) != null &&
                        (GetAgent(x + 1, y).type == BONE || GetAgent(x + 1, y).type == LINING)) {

                    PERF.Set(x, y, 0.1 * PERF_DiffCoef);
                }

                if (GetAgent(x, y) != null &&
                        (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {

                    PERF.Set(x, y, 0.1 * PERF_DiffCoef);

                } else if (y != yDim - 1 &&
                        GetAgent(x, y + 1) != null &&
                        (GetAgent(x, y + 1).type == BONE || GetAgent(x, y + 1).type == LINING)) {

                    PERF.Set(x, y, 0.1 * PERF_DiffCoef);
                }
            }
        }

// PERF Diffusion
        PERF.DiffusionADI(PERF_DiffCoef);

// Natural Decay of PERF
        PERF.MulAll(PERF_decayRate);
        PERF.Update();

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
    public void DrawPERF(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else
                    vis.SetPix(x, y, HeatMapRGB(PERF.Get(x, y)));

            }
        }
    }

    public void RecordOut(FileIO writeHere,int time){
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int[] cts = new int[9];

        for (simpleBoneCell c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            } else if(c.type==LINING){
                //ct_LINING++;
                cts[1]++;
            } else if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
                cts[2]++;
            } else if(c.type==MM && c.RESISTANT) {
                cts[3]++;
            } else if(c.type==MM && c.bcmaLoss) {
                cts[4]++;
            } else if(c.type==activeTcell) {
                cts[5]++;
            } else if(c.type==EXHT_CELL) {
                cts[6]++;
            } else if(c.type == supressorTcell){
                cts[7]++;
            }else if(c.type == naiveTcell){
                cts[8]++;
            }

        }
        //population of one timestep per line
        writeHere.Write(time+",");
        writeHere.WriteDelimit(cts,",");
        writeHere.Write("\n");
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
           if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
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
            }else if(c.type == BONE){
                cellType = "BONE";
            }else if(c.type == LINING){
                cellType = "LINING";
            }
            writeHere.Write(time + "," + simID + "," + position +  "," + cellType +"\n");
        }
    }

    public double[] CellCounts(){
        double[] cts = new double[9];

        for (simpleBoneCell c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            }else if(c.type==LINING){
                //ct_LINING++;
                cts[1]++;
            } else if(c.type==MM && !c.RESISTANT && !c.bcmaLoss){
                cts[2]++;
            } else if(c.type==MM && c.RESISTANT){
                cts[3]++;
            } else if(c.type==MM && c.bcmaLoss) {
                cts[4]++;
            } else if(c.type == activeTcell){
                cts[5]++;
            } else if(c.type == EXHT_CELL){
                cts[6]++;
            } else if(c.type == supressorTcell){
                cts[7]++;
            }else if(c.type == naiveTcell){
                cts[8]++;
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

    double tcellAge=0;
    double lifespan=0;
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
    boolean TGFB_on = false;
    boolean PERF_on;
    int cloneID;
    double p_kill = ProbScale(.9, TIMESTEP_AGENT);
    double extp_kill = ProbScale(.09, TIMESTEP_AGENT);
    int cellID;
    int timeSinceDivision;
    static final int[] MOORE_HOOD = MooreHood(true);
    static final int[] MOORE_HOOD_NO_CENTER = MooreHood(false);
    static final int[] VN_HOOD = VonNeumannHood(false);

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

    public int noisyseekCXCL9() {
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

    public int seekTGFB() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] TGFB_levels = new double[9]; // stores TGFB levels at the 8 directions and the center

        // Extracting TGFB levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            TGFB_levels[i] = G.TGFB.Get(G.tmoveHood[i]);
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
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[1] - TGFB_levels[2]);
                        break;
                    case 2: // left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[1] - TGFB_levels[2]);
                        break;
                    case 3: // up
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[3] - TGFB_levels[4]);
                        break;
                    case 4: // down
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[3] - TGFB_levels[4]);
                        break;
                    case 5: // top-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[5] - TGFB_levels[6]);
                        break;
                    case 6: // top-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[5] - TGFB_levels[6]);
                        break;
                    case 7: // bottom-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[7] - TGFB_levels[8]);
                        break;
                    case 8: // bottom-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[7] - TGFB_levels[8]);
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

    public int seekPerf() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] CXCL9_levels = new double[9]; // stores CXCL9 levels at the 8 directions and the center

        // Extracting CXCL9 levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.PERF.Get(G.tmoveHood[i]);
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

    public void CellStep(int time, double [] Cell_Counts, int simID) {

        ///////////
        //MYELOMA//
        ///////////
        if (type == MM) {

            if (time ==0){

                this.simulationID = simID;
            }

            double rn_BirthDeath = G.rn.Double();
            double pdiv;
            double pdeath = G.MM_DEATH;
            //int[] Hood = MooreHood(true); // For division and movement
            int options = MapOccupiedHood(MOORE_HOOD);
            double x = ProbScale(5.5e-4, TIMESTEP_AGENT);
            for (int j = 0; j < options; j++) {
                simpleBoneCell neighbor = G.GetAgent(MOORE_HOOD[j]);
                if (neighbor != null && (neighbor.type == BONE || neighbor.type == LINING)) {
                    if(G.rn.Double() < x){
                        TGFB_on=true;
                        neighbor.Dispose();
                        break;
                    }
                }
            }
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

                        if (G.rn.Double() < G.antigenLoss && this.mhcLoss==false){
                            double newExpression = 0;
                            child.mhcIExpression = newExpression;
                            child.mhcLoss = true;
                        }
                        else{child.mhcIExpression = this.mhcIExpression;
                        }

                        // Antigen loss mutation
                        if (G.rn.Double() < G.antigenLoss && this.bcmaLoss==false) {
                            double newExpression = 0;
                            if (child != null) {
                                child.bcmaExpression = newExpression;
                                child.bcmaLoss = true;
                            }
                        } else {
                            child.bcmaExpression = this.bcmaExpression;
                            child.mhcLoss = this.mhcLoss;
                        }
                    }
                }
            }
        }

        if (type == naiveTcell) {

            this.myeloma_bound = false;
            int[] hood = MooreHood(true);
            // Increment T cell age every 24 hours
            if (timeStep % (24 * 60 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                this.tcellAge += 1;
            }
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(hood);
            if (emptyNeighbors > 0 && G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) && this.timeSinceDivision >= 12) {
                this.timeSinceDivision = 0;
                int chosenCell = hood[G.rn.Int(emptyNeighbors)];
                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = naiveTcell;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                child.pd_l1 = child.pd_1 = G.boundedGaussian(2000, 1, 10, 2000);
            }
            for (int run = 0; run < 3; run++) {

                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots

                // Killing logic
                for (int j = 0; j < options; j++) {
                    if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM && G.GetAgent(movdivHood[j]).bcmaExpression > 0) {
                        if (G.rn.Double() < p_kill){
                            G.GetAgent(movdivHood[j]).Tcell_Kill();
                            this.pd_l1 = G.boundedGaussian(10, 1, 1, 20);
                            this.type = activeTcell;}
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
            if(G.TGFB.Get(Isq()) >= G.TGFBthresh){
                this.pd_1+=1;
            }
            // Exhaustion check
            if (this.pd_1 > pd_l1) {
                this.pd_1 = 0;
                this.pd_l1 = G.boundedGaussian(4, 1, 1, 22);
                this.type = EXHT_CELL;
                return;
            }

            // -----------------------------
            // Proliferation
            // -----------------------------
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(hood);
            if (emptyNeighbors > 0 && G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) && this.timeSinceDivision >= 24) {
                int chosenCell = hood[G.rn.Int(emptyNeighbors)];
                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = activeTcell;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(1, 1, 1, 2);
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
                    double kill_prob = 0;
                    if (target.mhcIExpression > 0) {
                        kill_prob = p_kill + (p_kill * 0.1);
                    }
                    if (G.rn.Double() < kill_prob) {
                        target.Tcell_Kill();
                        PERF_on = true;
                        this.pd_1 += 1;
                        break;
                    }else{
                        PERF_on = false;
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
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(hood);
            if (emptyNeighbors > 0 && G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) && this.timeSinceDivision >= 12) {
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
                            this.pd_1 += 1;
                            PERF_on = true;
                        }
                        else{
                            PERF_on = false;
                        }
                    }

                }

                // -----------------------------
                // Movement (dysfunctional chemotaxis)
                // -----------------------------
                int moveTo = noisyseekCXCL9();
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
            if (this.pd_1 >= this.pd_l1) {
                this.Dispose();
                return;
            }
            // Proliferation
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(hood);
            if (emptyNeighbors > 0 && G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) && this.timeSinceDivision >= 24) {
                int chosenCell = hood[G.rn.Int(emptyNeighbors)];
                simpleBoneCell child = G.NewAgentSQ(chosenCell);
                child.type = supressorTcell;
                child.tcellAge = 0;
                child.lifespan = G.boundedGaussian(30, 1, 30, 34);
                child.pd_l1 = G.boundedGaussian(20, 1, 10, 22);
            }
            // -----------------------------
            // Check neighbors for T cells
            // -----------------------------
            int occupied = MapOccupiedHood(hood);
            simpleBoneCell target = null;
            for (int j = 0; j < occupied; j++) {

                simpleBoneCell neighbor = G.GetAgent(hood[j]);

                if (neighbor.type == activeTcell ||
                        neighbor.type == EXHT_CELL ||
                        neighbor.type == naiveTcell) {

                    target = neighbor;
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
            int moveTo = seekPerf();
            if (moveTo >= 0 && G.GetAgent(moveTo) == null) {
                MoveSQ(moveTo);
            }
        }
    }
}
