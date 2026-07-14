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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     GRID CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public class simpleBoneGrid extends AgentGrid2D<simpleBoneCell> implements SerializableModel {

    ///////////////
    //GRID FIELDS//
    ///////////////

    // ----------------------------------------------------------------------------------------
    // DOSING SCHEDULE
    // A "Drug" is applied one per calendar week (see runSimulation). NONE represents a
    // treatment-free interval (drug holiday). Replaces the old pair of booleans
    // (BCMA_TCE / GPRC5D_TCE) that were previously toggled directly in scattered places, and
    // the old combo_tce boolean + int[] schedule (1/2/0 magic numbers) that only handled
    // switching between the two drugs and had no clean way to express "off".
    // ----------------------------------------------------------------------------------------
    public enum Drug { NONE, BCMA, GPRC5D }

    // Optional helper for building a schedule out of named blocks instead of writing out every
    // week by hand, e.g. "2 weeks BCMA, 4 week holiday, 2 weeks GPRC5D":
    //
    //   List<DoseBlock> plan = List.of(
    //       new DoseBlock(Drug.BCMA, 2),
    //       new DoseBlock(Drug.NONE, 4),
    //       new DoseBlock(Drug.GPRC5D, 2)
    //   );
    //   Drug[] schedule = expand(plan);
    //
    public record DoseBlock(Drug drug, int weeks) {}

    public static Drug[] expand(List<DoseBlock> plan) {
        List<Drug> out = new ArrayList<>();
        for (DoseBlock b : plan) {
            for (int i = 0; i < b.weeks(); i++) out.add(b.drug());
        }
        return out.toArray(new Drug[0]);
    }

    // ----------------------------------------------------------------------------------------
    // CONFIGURATION FLAGS
    // These are set once (by hand, before compiling) and only ever READ during a run, never
    // written mid-simulation, so they are safe to leave as static/shared across threads.
    // ----------------------------------------------------------------------------------------
    public static boolean TCELL = true;
    public static boolean BIWEEKLY = false;
    public static boolean ADAPTIVE = false;
    public static boolean TREATMENT_ON = false; //this is to control TGFBi

    //CLUSTER
    public static boolean PARAM_SWEEP = false; //use when importing parameters to loop through
    public static boolean runPar = true;
    public static boolean HEADLESS =true; //use true with cluster
    public static boolean GIFSAVE = false;
    public static boolean LOCAL = true; // use false with cluster
    public static double numSteps = 1.0 * 365.0 * 24.0 * 60.0; // years the model will run
    public static int numSims = 10; //Number of Simulations

    // ----------------------------------------------------------------------------------------
    // RUNTIME TREATMENT STATE
    // IMPORTANT: these used to be `static`, which is a bug once simulations run concurrently on
    // a thread pool (every thread would read/write the SAME field, corrupting each other's
    // dosing schedules). They are now per-grid instance fields, so every simulation gets its
    // own independent copy.
    // ----------------------------------------------------------------------------------------
    public boolean TCE = true;
    public boolean BCMA_TCE = true;
    public boolean GPRC5D_TCE = false; // was GRPC5_TCE (also fixes the GPRC5D <-> GRPC5D typo)

    // Unique, ever-increasing clone ID counter for this simulation's lineage tracking.
    // Previously `cloneID` was incremented on the parent cell itself at division time
    // (`child.cloneID = this.cloneID += 1`), which mutated the parent's own ID and did not
    // produce stable, unique per-cell IDs. Now every new cell gets a fresh unique ID from this
    // counter, and `parentID` records the ID of the cell it divided from.
    public int nextCloneID = 1;

    // Set once per simulation (see runSimulation) and used for all output records, instead of
    // relying on each individual cell to record its own simulationID (which was previously only
    // ever set on MM founder cells at time 0, leaving every other cell/row defaulted to 0).
    public int simulationID;

    public final static int BONE = RGB256(255, 255, 250),
            LINING = RGB256(64, 106, 151),
            MM = RGB256(0, 128, 0),
            activeTcell = RGB256(17, 150, 150),
            EXHT_CELL = RGB256(200, 50, 250),
            supressorTcell = RGB256(255, 165, 0),
            naiveTcell = RGB256(100, 100, 100),
            extnaiveTcell = RGB256(100, 115, 100);

    // ----------------------------------------------------------------------------------------
    // Cell_Counts / cts array layout — named indices so the array isn't addressed by magic
    // numbers in three different places (CellCounts, RecordOut, the ADAPTIVE-therapy check).
    // ----------------------------------------------------------------------------------------
    public static final int IDX_BONE = 0;
    public static final int IDX_LINING = 1;
    public static final int IDX_MM_SENSITIVE = 2;
    public static final int IDX_MM_RESISTANT = 3;
    public static final int IDX_MM_BCMA_LOSS = 4;
    public static final int IDX_TCELL_ACTIVE = 5;
    public static final int IDX_TCELL_EXHAUSTED = 6;
    public static final int IDX_TCELL_TREG = 7;
    public static final int IDX_TCELL_NAIVE = 8;
    public static final int IDX_TCELL_EXT_NAIVE = 9;
    public static final int IDX_MM_MHC_LOSS = 10;
    public static final int IDX_MM_GPRC5D_LOSS = 11;
    public static final int IDX_MM_DOUBLE_LOSS = 12;
    public static final int CELL_COUNT_ARRAY_SIZE = 13;

    //SETUP
    static double MinToHour = 60.0;
    public final static double SPACESTEP = 10.0;//um
    public static double TIMESTEP_AGENT = 6.0 / MinToHour; //0.1;//hr; //6.0 min; 6.0/60.0 hour
    public final static double N_TIMESTEP_PDE = 60.0 * (MinToHour * TIMESTEP_AGENT);//360.0; //Number of diffusion timesteps; 1 dts = 1 sec; 360 dts = 1 ts = 6 min
    // How many agent timesteps make up one day — computed once instead of recomputed inline
    // at every call site (main loop daily checks, GIF-frame cadence, per-cell aging).
    public static final double TIMESTEPS_PER_DAY = 24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT);

    //DiffCoef MUST <0.25 for FTCS scheme!

    //CHEMOTAXIS
    double Tcell_DiffCoef = 0.01 * 3.0 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP);
    double Tcell_TaxisCoeff = 5.0e10 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP);

    //CELL PARAMETERS
    public int TURNOVER_TIME = (int) (2102400.0 / (MinToHour * TIMESTEP_AGENT)); //2102400 min = 350400 ts = 4 years
    public double MM_DEATH = 1.0 / 11000 * (MinToHour);
    public double pmutate = 0.0;
    public double antigenLoss = Math.pow(10, -4);
    public double T_CELL_DIV_RATE = 1.0 / 1440 * (MinToHour); // T_CELL DIVISION RATE
    double CXCL9_productionRate = (2.04e-9 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE) / 8;
    double CXCL9_decayRate = -.2 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE;
    double CXCL9_DiffCoef = 2700.0 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP * N_TIMESTEP_PDE);
    double maxCXCL9 = (2.07e-9) / 8;

    double PERF_productionRate = 2.04e-9 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE;
    static double PERF_basalRate = 2.04e-11 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE;
    static double PERF_decayRate = -0.35 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE; // BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    double PERF_DiffCoef = 780.0 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP * N_TIMESTEP_PDE);
    double Extra_PERFtime = 4320.0 / (MinToHour * TIMESTEP_AGENT); // 4320 min = 720 ts = 3 days
    static double maxPERF = 8.7e-10;

    //Cytokines
    double TGFB_productionRate = 2.04e-9 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE;
    static double TGFB_basalRate = 2.04e-11 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE;
    static double TGFB_decayRate = -0.35 * (MinToHour * TIMESTEP_AGENT) / N_TIMESTEP_PDE; // BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    double TGFB_DiffCoef = 780.0 * (MinToHour * TIMESTEP_AGENT) / (SPACESTEP * SPACESTEP * N_TIMESTEP_PDE);
    double Extra_TGFBtime = 4320.0 / (MinToHour * TIMESTEP_AGENT); // 4320 min = 720 ts = 3 days
    static double maxTGFB = 8.7e-10;

    double Ts = TGFB_basalRate / (Math.abs(TGFB_decayRate) * maxTGFB); //Basal TGFB
    double TGFBthresh = (1.05) * TGFB_basalRate / (Math.abs(TGFB_decayRate) * maxTGFB); //was 0.01; now 5% increase from basal

    //MODEL TESTS
    public double MarrowArea;
    double convert_to_days = (MinToHour * TIMESTEP_AGENT) / (60.0 * 24.0); //1 ts = 6 min = 1/240 day
    int count_BA = 0;
    int init_BA = 0;
    int Nts = (int) ((numSteps) / (MinToHour * TIMESTEP_AGENT));

    public Rand rn;
    public PDEGrid2D CXCL9;
    public PDEGrid2D TGFB;
    public PDEGrid2D PERF;

    public int[] tmoveHood = MooreHood(true);
    public ArrayList<Integer> InitBoneList = new ArrayList<>();
    public ArrayList<simpleBoneCell> AllBoneList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<simpleBoneCell> LiningList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs

    FileIO out;
    FileIO clones;
    FileIO bones;
    FileIO locations;
    FileIO InitialBone;
    FileIO tgfbLocs;

    //This is important for serializable model
    @Override
    public void SetupConstructors() {
        this._PassAgentConstructor(simpleBoneCell.class);
    }

    ////////////////////
    //GRID CONSTRUCTOR//
    ////////////////////

    public simpleBoneGrid(int xDim, int yDim, Rand rn, String Bone_FileName) {
        super(xDim, yDim, simpleBoneCell.class, true, true);
        this.rn = rn;

        //Create 2D PDE Grid for RANKL/CXCL9/TGFB/PERF fields and boundary condition
        CXCL9 = new PDEGrid2D(xDim, yDim, true, true);
        TGFB = new PDEGrid2D(xDim, yDim, true, true);
        PERF = new PDEGrid2D(xDim, yDim, true, true);
        InitialBone = new FileIO(Bone_FileName, "r");
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

    public static void runSimulation(final int simID, final int prow, final String baseFolder, final ArrayList<String> param_list, final Drug[] schedule) {

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
        UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
        UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 2);
        UIGrid TGFB_vis = new UIGrid(xDim, yDim, 2);
        UIGrid PERF_vis = new UIGrid(xDim, yDim, 2);
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
        // GIF makers — only stood up when they'll actually be used, so headless/no-GIF runs
        // don't pay for GIF file setup they never touch.
        // -------------------------
        GifMaker gm_Cell_vis = null;
        GifMaker gm_CXCL9_vis = null;
        GifMaker gm_TGFB_vis = null;
        GifMaker gm_PERF_vis = null;
        if (!HEADLESS && GIFSAVE) {
            gm_Cell_vis = new GifMaker(subfolder.concat("/CellVid.gif"), 0, true);
            gm_CXCL9_vis = new GifMaker(subfolder.concat("/CXCL9.gif"), 0, true);
            gm_TGFB_vis = new GifMaker(subfolder.concat("/TGFB.gif"), 0, true);
            gm_PERF_vis = new GifMaker(subfolder.concat("/PERF.gif"), 0, true);
        }

        // -------------------------
        // Bone file selection and grid creation
        // -------------------------
        String Bone_Filename;
        if (LOCAL) {
            Bone_Filename = "/Users/80024703/Desktop/code/Bone/BAout_2020May5_Sim14.csv";
        } else {
            Bone_Filename = "Bone/BAout_2020May5_Sim14.csv";
        }

        // Create grid and RNG (deterministic seed per sim/param-row combination)
        Rand rng = new Rand(simID + 1000L * prow);
        simpleBoneGrid g = new simpleBoneGrid(xDim, yDim, rng, Bone_Filename);
        g.simulationID = simID;

        // Set parameters from param_list if param sweep active
        if (PARAM_SWEEP && param_list != null) {
            g.SetParams(prow, param_list);
        }

        // Record Output (create PopOut.csv in this sim folder)
        g.newFileIO(subfolder, "w");

        // Initialize model
        g.InitBone();

        // -------------------------
        // main time loop
        // -------------------------
        double exhaustedFraction = 0.25;
        int totalTcell = 250;

        for (int i = 0; i < g.Nts; i++) {
            int daysPassedCalc = (int) Math.floor(i / TIMESTEPS_PER_DAY);
            double[] Cell_Counts = g.CellCounts();

            // ----------------------------------------------------------------------------
            // DOSING SCHEDULE
            // One Drug (or NONE) per calendar week, taken from `schedule` and cycled once the
            // run outlasts it. Replaces the old BIWEEKLY-toggle / combo_tce blocks: a schedule
            // of all-BCMA weeks reproduces old constant-BCMA behavior, and inserting Drug.NONE
            // entries gives treatment-free intervals for free. BIWEEKLY/ADAPTIVE below can
            // still override this if you need those specific behaviors for a run.
            // ----------------------------------------------------------------------------
            if (schedule != null && schedule.length > 0) {
                int week = daysPassedCalc / 7;
                Drug currentDrug = schedule[week % schedule.length];

                g.BCMA_TCE = (currentDrug == Drug.BCMA);
                g.GPRC5D_TCE = (currentDrug == Drug.GPRC5D);
                g.TCE = (currentDrug != Drug.NONE);
            }

            // Weekly therapy toggle (overrides the schedule above when enabled)
            if (BIWEEKLY) {
                if (daysPassedCalc % 7 == 0 && i % (int) TIMESTEPS_PER_DAY == 0) {
                    g.TCE = true;
                } else {
                    g.TCE = false;
                }
            }

            if (ADAPTIVE && (Cell_Counts[IDX_MM_SENSITIVE] + Cell_Counts[IDX_MM_RESISTANT] + Cell_Counts[IDX_MM_BCMA_LOSS]
                    + Cell_Counts[IDX_MM_MHC_LOSS] + Cell_Counts[IDX_MM_GPRC5D_LOSS] + Cell_Counts[IDX_MM_DOUBLE_LOSS]) <= 3000) {
                g.TCE = false;
            } else if (ADAPTIVE) {
                g.TCE = true;
            }

            if (!HEADLESS) {
                win.TickPause(10);
            }

            if (TCELL) {
                if (i == 0) {
                    seedTcells(g, totalTcell, exhaustedFraction, false);
                } else if ((Cell_Counts[IDX_TCELL_ACTIVE] + Cell_Counts[IDX_TCELL_EXHAUSTED]
                        + Cell_Counts[IDX_TCELL_NAIVE] + Cell_Counts[IDX_TCELL_EXT_NAIVE]) < totalTcell && g.TCE) {
                    int replenishCount = (BIWEEKLY || ADAPTIVE|| schedule != null) ? 250 : 10;
                    seedTcells(g, replenishCount, exhaustedFraction, true);
                }
            }

            // Model step and drawing
            g.ModelStep(i, Cell_Counts, simID);

            if (!HEADLESS) {
                g.Draw(Cell_vis, days, i, simID);
                g.DrawCXCL9(CXCL9_vis);
                g.DrawTGFB(TGFB_vis);
                g.DrawPERF(PERF_vis);
                if (GIFSAVE && i % (int) TIMESTEPS_PER_DAY == 0) {
                    gm_Cell_vis.AddFrame(Cell_vis);
                    gm_CXCL9_vis.AddFrame(CXCL9_vis);
                    gm_TGFB_vis.AddFrame(TGFB_vis);
                    gm_PERF_vis.AddFrame(PERF_vis);
                }
            }

            // daily recording
            if (i % (int) TIMESTEPS_PER_DAY == 0) {
                g.RecordOut(g.out, i);
                g.RecordClones(g.clones, i);
                g.RecordBones(g.bones, i);
                g.RecordLocs(g.locations, i);
                g.RecordTGFBLocs(g.tgfbLocs, i);
            }
        } // end time loop

        // -------------------------
        // cleanup: close files, gifs, window
        // -------------------------
        g.closeFileIO();
        if (gm_Cell_vis != null) {
            gm_Cell_vis.Close();
            gm_CXCL9_vis.Close();
            gm_TGFB_vis.Close();
            gm_PERF_vis.Close();
        }

        if (!HEADLESS) {
            win.Close();
        }
    }

    // ----------------------------------------------------------------------------------------
    // T-cell seeding helpers
    //
    // Consolidates what used to be four separate, near-duplicated ~20-line blocks (initial
    // random+boundary placement, and replenishment random+boundary placement) into one
    // parameterized helper. Also fixes a bug in the original: the random-site attempt and the
    // boundary-fallback attempt both ran unconditionally within the same loop iteration, so a
    // single "attempt" could place up to two T cells and overshoot the intended population by
    // as much as 2x. Now exactly one placement is attempted (and counted) per loop iteration.
    // Also replaces the unseeded java.util.Random used for placement with the grid's own seeded
    // `rn`, so runs are actually reproducible given a fixed seed.
    // ----------------------------------------------------------------------------------------
    private static boolean placeTcellAt(simpleBoneGrid g, int siteIndex, double exhaustedFraction) {
        if (g.GetAgent(siteIndex) != null) {
            return false;
        }
        simpleBoneCell c = g.NewAgentSQ(siteIndex);
        c.lifespan = g.boundedGaussian(30, 1, 30, 34);
        if (g.rn.Double() < .1) {
            c.type = supressorTcell;
            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
        } else if (g.rn.Double() < exhaustedFraction) {
            c.type = EXHT_CELL;
            c.pd_l1 = g.boundedGaussian(10, 1, 1, 20);
        } else {
            c.type = activeTcell;
            c.pd_l1 = g.boundedGaussian(20, 1, 1, 40);
        }
        return true;
    }

    private static void seedTcells(simpleBoneGrid g, int countToPlace, double exhaustedFraction, boolean boundaryOnly) {
        int placed = 0;
        int attempts = 0;
        int maxAttempts = 100000;
        int[] boundaries = g.BoundaryIs();

        while (placed < countToPlace && attempts < maxAttempts) {
            attempts++;
            int siteIndex;
            if (!boundaryOnly && placeTcellAt(g, g.rn.Int(g.length), exhaustedFraction)) {
                placed++;
                continue;
            }
            siteIndex = boundaries[g.rn.Int(boundaries.length)];
            if (placeTcellAt(g, siteIndex, exhaustedFraction)) {
                placed++;
            }
        }
    }

    //sample from a bounded distribution
    public double boundedGaussian(double mean, double dev, double min, double max) {
        double gauss = rn.Gaussian(0, 1);
        double val = dev * gauss + mean;
        while (val > max || val < min) {
            gauss = rn.Gaussian(0, 1);
            val = dev * gauss + mean;
        }
        return val;
    }

    public void newFileIO(String projPath, String mode) {

        out = new FileIO(projPath + "PopOut.csv", mode);
        clones = new FileIO(projPath + "clones.csv", mode);
        bones = new FileIO(projPath + "boneOut.csv", mode);
        locations = new FileIO(projPath + "cellLocs.csv", mode);
        tgfbLocs = new FileIO(projPath + "tgfbLoc.csv", mode);

        if (mode.equals("w")) {
            out.Write("Timestep" + ","
                    + "BONE" + ","
                    + "LINING" + ","
                    + "S_MM" + ","
                    + "R_MM" + ","
                    + "AL_MM" + ","
                    + "TCell" + ","
                    + "ExtTcell" + ","
                    + "T-reg" + ","
                    + "Naive Tcell" + ","
                    + "Ext Naive Tcell" + ","
                    + "ML_MM" + ","
                    + "GL_MM" + ","
                    + "AL_GL_MM" + ","
                    + "\n");
            clones.Write("Timestep" + "," + "parentID" + "," + "cloneID" + "," + "MHCI" + "," + "BCMA" + "," + "GPRC5D" + "\n");
            bones.Write("Timestep" + "," + "SimID" + "," + "Position" + "\n");
            locations.Write("Timestep" + "," + "SimID" + "," + "Position" + "," + "Type" + "\n");
            tgfbLocs.Write("TimeStep" + "," + "Position" + "," + "Concentration" + "\n");
        }
    }

    public void closeFileIO() {
        out.Close();
        clones.Close();
        bones.Close();
        locations.Close();
        tgfbLocs.Close();
    }

    public void SetParams(int prow, ArrayList<String> param_list) {
        String[] split_param_list = param_list.get(prow).split(",");
        pmutate = Double.parseDouble(split_param_list[2]);
    }

    public void InitBone() {
        int initMyeloma = 3000;

        //  FOR IRREGULAR BONE
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
            if (GetAgent(Integer.parseInt(split_input_data[index])).MarrowInHood()) {
                GetAgent(Integer.parseInt(split_input_data[index])).type = LINING;
                GetAgent(Integer.parseInt(split_input_data[index])).liningAge = TURNOVER_TIME;
                LiningList.add(GetAgent(Integer.parseInt(split_input_data[index])));
            }
        }
        init_BA = InitBoneList.size();
        MarrowArea = (xDim * yDim) - init_BA;
        int myelomaCellsToPlace = initMyeloma;
        int placedMyelomaCells = 0;
        int boneProximityDistance = 10; // Maximum initial distance from bone

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
                        double distance = Math.sqrt(Math.pow(xi - x, 2) + Math.pow(yi - y, 2));
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

            // Place myeloma cell as a new founder clone
            simpleBoneCell mm = NewAgentSQ(x, y);
            mm.type = MM;
            mm.bcmaExpression = 1;
            mm.gprc5dExpression = 1;
            mm.mhcIExpression = 1;
            mm.cloneID = nextCloneID++;
            mm.parentID = 0; // 0 marks a founder cell with no parent

            placedMyelomaCells++;
        }
    }

    public void ModelStep(int time, double[] Cell_Counts, int simID) {

        /////////////////////////////////////////////////
        //STEP 1: REACTION-DIFFUSION FOR CXCL9/TGFB/PERF//
        /////////////////////////////////////////////////

        for (int x = 0; x < CXCL9.xDim; x++) {
            for (int y = 0; y < CXCL9.yDim; y++) {
                simpleBoneCell here = GetAgent(x, y);
                simpleBoneCell rightNeighbor = (x != xDim - 1) ? GetAgent(x + 1, y) : null;
                simpleBoneCell upNeighbor = (y != yDim - 1) ? GetAgent(x, y + 1) : null;

                if (here != null && here.type == MM && here.bcmaExpression > 0) {
                    CXCL9.Add(x, y, CXCL9_productionRate / maxCXCL9);
                }

                boolean nearBoneOrLining =
                        (here != null && (here.type == BONE || here.type == LINING)) ||
                                (rightNeighbor != null && (rightNeighbor.type == BONE || rightNeighbor.type == LINING)) ||
                                (upNeighbor != null && (upNeighbor.type == BONE || upNeighbor.type == LINING));

                if (nearBoneOrLining) {
                    CXCL9.Set(x, y, 0.1 * CXCL9_DiffCoef);
                }
            }
        }

        CXCL9.DiffusionADI(CXCL9_DiffCoef);
        CXCL9.MulAll(CXCL9_decayRate);
        CXCL9.Update();

        for (int x = 0; x < TGFB.xDim; x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if (!TREATMENT_ON) {
                    simpleBoneCell here = GetAgent(x, y);
                    TGFB.Add(x, y, TGFB_basalRate / maxTGFB);
                    if (here != null && here.type == MM && here.TGFB_on) {
                        TGFB.Add(x, y, TGFB_productionRate / maxTGFB);
                    }
                }
                // when TREATMENT_ON, TGFB production is fully suppressed (nothing to add)
            }
        }

        TGFB.DiffusionADI(TGFB_DiffCoef);
        TGFB.MulAll(TGFB_decayRate);
        TGFB.Update();

        for (int x = 0; x < PERF.xDim; x++) {
            for (int y = 0; y < PERF.yDim; y++) {
                simpleBoneCell here = GetAgent(x, y);
                if (here != null && here.PERF_on) {
                    PERF.Add(x, y, PERF_productionRate / maxPERF);
                }
            }
        }

        PERF.DiffusionADI(PERF_DiffCoef);
        PERF.MulAll(PERF_decayRate);
        PERF.Update();

        /////////////////////////////////
        //STEP 2: ITERATE THROUGH CELLS//
        /////////////////////////////////

        CleanShuffle(rn);
        for (simpleBoneCell c : this) {
            c.CellStep(time, Cell_Counts, simID);
        }
    }

    ////////////////////////////
    //Full Domain (No zoom-in)//
    ////////////////////////////
    public void Draw(UIGrid vis, UILabel days, int i, int simID) {
        days.SetText("days: " + i * convert_to_days);
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x, y);
                if (drawMe != null && drawMe.type == MM && drawMe.RESISTANT) {
                    vis.SetPix(x, y, BLACK);
                } else if (drawMe != null && drawMe.type == MM && drawMe.bcmaLoss) {
                    vis.SetPix(x, y, BLACK);
                } else if (drawMe != null) {
                    vis.SetPix(x, y, drawMe.type);
                } else {
                    vis.SetPix(x, y, RGB256(240, 220, 220)); //MARROW=LIGHT PINK
                }
            }
        }
        vis.SetString("Day: " + (int) (i * convert_to_days), 1, yDim - 1, BLACK, RGB256(240, 220, 220));
    }

    public void DrawCXCL9(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x, y);
                if (drawMe != null && drawMe.type == LINING) {
                    vis.SetPix(x, y, BONE);
                } else {
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));
                }
            }
        }
    }

    public void DrawTGFB(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x, y);
                if (drawMe != null && drawMe.type == LINING) {
                    vis.SetPix(x, y, BONE);
                } else {
                    vis.SetPix(x, y, HeatMapRGB(TGFB.Get(x, y)));
                }
            }
        }
    }

    public void DrawPERF(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                simpleBoneCell drawMe = GetAgent(x, y);
                if (drawMe != null && drawMe.type == LINING) {
                    vis.SetPix(x, y, BONE);
                } else {
                    vis.SetPix(x, y, HeatMapRGB(PERF.Get(x, y)));
                }
            }
        }
    }

    public void RecordOut(FileIO writeHere, int time) {
        int[] cts = new int[CELL_COUNT_ARRAY_SIZE];

        for (simpleBoneCell c : this) {
            if (c.type == BONE) {
                cts[IDX_BONE]++;
            } else if (c.type == LINING) {
                cts[IDX_LINING]++;
            } else if (c.type == activeTcell) {
                cts[IDX_TCELL_ACTIVE]++;
            } else if (c.type == EXHT_CELL) {
                cts[IDX_TCELL_EXHAUSTED]++;
            } else if (c.type == supressorTcell) {
                cts[IDX_TCELL_TREG]++;
            } else if (c.type == naiveTcell) {
                cts[IDX_TCELL_NAIVE]++;
            } else if (c.type == extnaiveTcell) {
                cts[IDX_TCELL_EXT_NAIVE]++;
            } else if (c.type == MM) {
                if (c.gprc5dLoss && c.bcmaLoss) {
                    cts[IDX_MM_DOUBLE_LOSS]++;
                } else if (c.gprc5dLoss) {
                    cts[IDX_MM_GPRC5D_LOSS]++;
                } else if (c.bcmaLoss) {
                    cts[IDX_MM_BCMA_LOSS]++;
                } else if (c.mhcLoss) {
                    cts[IDX_MM_MHC_LOSS]++;
                } else if (c.RESISTANT) {
                    cts[IDX_MM_RESISTANT]++;
                } else {
                    cts[IDX_MM_SENSITIVE]++;
                }
            }
        }
        //population of one timestep per line
        writeHere.Write(time + ",");
        writeHere.WriteDelimit(cts, ",");
        writeHere.Write("\n");
    }

    public void RecordClones(FileIO writeHere, int time) {
        for (simpleBoneCell c : this) {
            if (c.type == MM) {
                writeHere.Write(time + "," + c.parentID + "," + c.cloneID + "," + c.mhcIExpression + ","
                        + c.bcmaExpression + "," + c.gprc5dExpression + "\n");
            }
        }
    }

    public void RecordBones(FileIO writeHere, int time) {
        for (simpleBoneCell c : this) {
            if (c.type == BONE) {
                writeHere.Write(time + "," + simulationID + "," + c.Isq() + "\n");
            }
        }
    }

    public void RecordLocs(FileIO writeHere, int time) {
        for (simpleBoneCell c : this) {
            String cellType = "UNKNOWN";
            if (c.type == MM) {
                if (c.bcmaLoss && c.gprc5dLoss) {
                    cellType = "AL_GL_MM";
                } else if (c.bcmaLoss) {
                    cellType = "AL_MM";
                } else if (c.gprc5dLoss) {
                    cellType = "GL_MM";
                } else if (c.mhcLoss) {
                    cellType = "MHC_LOSS_MM";
                } else if (c.RESISTANT) {
                    cellType = "Resistant MM";
                } else {
                    cellType = "Sensitive MM";
                }
            } else if (c.type == activeTcell) {
                cellType = "Tcell";
            } else if (c.type == EXHT_CELL) {
                cellType = "Ext_Tcell";
            } else if (c.type == supressorTcell) {
                cellType = "Treg";
            } else if (c.type == BONE) {
                cellType = "BONE";
            } else if (c.type == LINING) {
                cellType = "LINING";
            }
            writeHere.Write(time + "," + simulationID + "," + c.Isq() + "," + cellType + "\n");
        }
    }

    public void RecordTGFBLocs(FileIO writeHere, int time) {
        for (int x = 0; x < TGFB.xDim; x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                double concentration = TGFB.Get(x, y);
                int position = x * TGFB.yDim + y;
                writeHere.Write(time + "," + position + "," + concentration + "\n");
            }
        }
    }

    public double[] CellCounts() {
        double[] cts = new double[CELL_COUNT_ARRAY_SIZE];
        for (simpleBoneCell c : this) {
            if (c.type == BONE) {
                cts[IDX_BONE]++;
            } else if (c.type == LINING) {
                cts[IDX_LINING]++;
            } else if (c.type == MM) {
                if (c.bcmaLoss && c.gprc5dLoss) {
                    cts[IDX_MM_DOUBLE_LOSS]++;
                } else if (c.bcmaLoss) {
                    cts[IDX_MM_BCMA_LOSS]++;
                } else if (c.gprc5dLoss) {
                    cts[IDX_MM_GPRC5D_LOSS]++;
                } else if (c.mhcLoss) {
                    cts[IDX_MM_MHC_LOSS]++;
                } else if (c.RESISTANT) {
                    cts[IDX_MM_RESISTANT]++;
                } else {
                    cts[IDX_MM_SENSITIVE]++;
                }
            } else if (c.type == activeTcell) {
                cts[IDX_TCELL_ACTIVE]++;
            } else if (c.type == EXHT_CELL) {
                cts[IDX_TCELL_EXHAUSTED]++;
            } else if (c.type == supressorTcell) {
                cts[IDX_TCELL_TREG]++;
            } else if (c.type == naiveTcell) {
                cts[IDX_TCELL_NAIVE]++;
            } else if (c.type == extnaiveTcell) {
                cts[IDX_TCELL_EXT_NAIVE]++;
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

        // ----------------------------------------------------------------------------------
        // DOSING SCHEDULE(S) TO TEST
        // Add/remove entries here to sweep multiple drug schedules against every parameter row
        // and every replicate. Each schedule gets its own labeled output subfolder.
        // Example schedules:
        //   - constant BCMA every week
        //   - alternating BCMA/GPRC5D weekly
        //   - 2 weeks on, 2 weeks off (drug holiday)
        // ----------------------------------------------------------------------------------
        Map<String, Drug[]> schedulesToTest = new LinkedHashMap<>();
//        schedulesToTest.put("bcma_only", new Drug[]{ Drug.BCMA });
        schedulesToTest.put("alternating_weekly", new Drug[]{ Drug.BCMA, Drug.GPRC5D });
//        schedulesToTest.put("alt_bcma_gprc5d_on1_off1", expand(List.of(
//                new DoseBlock(Drug.BCMA, 1),
//                new DoseBlock(Drug.NONE, 1),
//                new DoseBlock(Drug.GPRC5D, 1),
//                new DoseBlock(Drug.NONE, 1)
//
//        )));

        // Thread pool (limit threads if you want)
        int numThreads = Math.min(numSims, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        if (runPar) {
            for (Map.Entry<String, Drug[]> scheduleEntry : schedulesToTest.entrySet()) {
                final String scheduleName = scheduleEntry.getKey();
                final Drug[] schedule = scheduleEntry.getValue();
                final String scheduleFolder = fn + "/" + scheduleName;
                new File(scheduleFolder).mkdirs();

                for (int prow = 1; prow < param_list_size; prow++) {
                    for (int sim = 0; sim < numSims; sim++) {
                        final int simID = sim;
                        final int row = prow;
                        final ArrayList<String> paramsCopy = param_list; // may be null

                        executor.submit(() -> runSimulation(simID, row, scheduleFolder, paramsCopy, schedule));
                    }
                }
            }
            executor.shutdown();
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            System.out.println("All simulations finished.");
        } else {
            for (Map.Entry<String, Drug[]> scheduleEntry : schedulesToTest.entrySet()) {
                final String scheduleName = scheduleEntry.getKey();
                final Drug[] schedule = scheduleEntry.getValue();
                final String scheduleFolder = fn + "/" + scheduleName;
                new File(scheduleFolder).mkdirs();

                for (int prow = 1; prow < param_list_size; prow++) {
                    for (int sim = 0; sim < numSims; sim++) {
                        final int simID = sim;
                        final int row = prow;
                        final ArrayList<String> paramsCopy = param_list; // may be null
                        runSimulation(simID, row, scheduleFolder, paramsCopy, schedule);
                    }
                }
            }
            executor.shutdown();
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
    boolean gprc5dLoss = false; // was grpc5dLoss (GPRC5D <-> GRPC5D typo fix)
    boolean mhcLoss = false;
    double bcmaExpression = 1.0;
    double mhcIExpression = 1.0;
    double gprc5dExpression = 1.0; // was grpc5dExpression

    // Lineage tracking: cloneID is a unique ID assigned to this cell when it was created;
    // parentID is the cloneID of the cell it divided from (0 = founder cell, no parent).
    int parentID;
    int cloneID;

    //T cells
    double tcellAge = 0;
    double lifespan = 0;
    double pd_1 = 0;
    double pd_l1 = 0;
    boolean myeloma_bound = false;
    boolean TGFB_on = false;
    boolean PERF_on = false;

    double p_kill = ProbScale(.9, TIMESTEP_AGENT);
    double extp_kill = ProbScale(.3, TIMESTEP_AGENT);
    int timeSinceDivision;

    static final int[] MOORE_HOOD = MooreHood(true);
    static final int[] MOORE_HOOD_NO_CENTER = MooreHood(false);
    static final int[] VN_HOOD = VonNeumannHood(false);

    void Tcell_Kill() {
        if (type == MM) {
            this.Dispose();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CELL METHODS//////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    public int seekCXCL9() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] CXCL9_levels = new double[9];

        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.CXCL9.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>();

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]);

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

        double effectiveTaxis = G.Tcell_TaxisCoeff * (1.0 - exhaustionLevel);
        double noiseStrength = exhaustionLevel * G.Tcell_TaxisCoeff;
        double randomMoveProb = 0.2 * exhaustionLevel;

        if (G.rn.Double() < randomMoveProb) {
            return G.tmoveHood[G.rn.Int(neighbors)];
        }

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]);

                double P;
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

                double noisyGradient = gradientTerm + G.rn.Double() * noiseStrength;

                P = G.Tcell_DiffCoef + (effectiveTaxis * G.maxCXCL9 / 8.0) * noisyGradient;
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
        int neighbors = MapHood(G.tmoveHood);
        double[] TGFB_levels = new double[9];

        for (int i = 0; i < 9; i++) {
            TGFB_levels[i] = G.TGFB.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>();

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]);

                double P = 0;
                switch (i) {
                    case 1:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[1] - TGFB_levels[2]);
                        break;
                    case 2:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[1] - TGFB_levels[2]);
                        break;
                    case 3:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[3] - TGFB_levels[4]);
                        break;
                    case 4:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[3] - TGFB_levels[4]);
                        break;
                    case 5:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[5] - TGFB_levels[6]);
                        break;
                    case 6:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[5] - TGFB_levels[6]);
                        break;
                    case 7:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[7] - TGFB_levels[8]);
                        break;
                    case 8:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxTGFB) / 8 * (TGFB_levels[7] - TGFB_levels[8]);
                        break;
                }

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

        for (int i = 0; i < cProbArray.size(); i++) {
            if (rnum <= cProbArray.get(i) / ProbSum) {
                moveToIndex = emptyHood.get(i);
                break;
            }
        }

        return moveToIndex;
    }

    public int seekPerf() {
        int neighbors = MapHood(G.tmoveHood);
        double[] PERF_levels = new double[9];

        for (int i = 0; i < 9; i++) {
            PERF_levels[i] = G.PERF.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>();

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]);

                double P = 0;
                switch (i) {
                    case 1:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[1] - PERF_levels[2]);
                        break;
                    case 2:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[1] - PERF_levels[2]);
                        break;
                    case 3:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[3] - PERF_levels[4]);
                        break;
                    case 4:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[3] - PERF_levels[4]);
                        break;
                    case 5:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[5] - PERF_levels[6]);
                        break;
                    case 6:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[5] - PERF_levels[6]);
                        break;
                    case 7:
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[7] - PERF_levels[8]);
                        break;
                    case 8:
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (PERF_levels[7] - PERF_levels[8]);
                        break;
                }

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

        for (int i = 0; i < cProbArray.size(); i++) {
            if (rnum <= cProbArray.get(i) / ProbSum) {
                moveToIndex = emptyHood.get(i);
                break;
            }
        }

        return moveToIndex;
    }

    void Init() {
        if (type == BONE || type == LINING) {
            G.count_BA++;
        }
    }

    public boolean MarrowInHood() {
        int[] MarrowHood = VonNeumannHood(false); //4 neighbors
        int options = MapEmptyHood(MarrowHood);
        return options > 0;
    }

    double TGFBthresh = (1.05) * TGFB_basalRate / (Math.abs(TGFB_decayRate) * maxTGFB);

    public double Prob_MM_Divide(double TGFBval, double MaxDivRate, double halfmax) {
        double n = 2;
        return MaxDivRate * (1 / (1 + Math.pow((halfmax / TGFBval), n)));
    }

    public void CellStep(int time, double[] Cell_Counts, int simID) {

        ///////////
        //MYELOMA//
        ///////////
        if (type == MM) {

            double rn_BirthDeath = G.rn.Double();
            double pdeath = G.MM_DEATH;
            int options = MapOccupiedHood(MOORE_HOOD);
            double x = ProbScale(5.5e-4, TIMESTEP_AGENT);
            for (int j = 0; j < options; j++) {
                simpleBoneCell neighbor = G.GetAgent(MOORE_HOOD[j]);
                if (neighbor != null && (neighbor.type == BONE || neighbor.type == LINING) && G.rn.Double() < x) {
                    TGFB_on = true;
                    neighbor.Dispose();
                    break;
                }
            }
            double scaleFactor = 0.7 + (0.3 * this.bcmaExpression);
            // NOTE: previously this rate was run through ProbScale() twice — once here to turn
            // the per-day division rate into an hourly-equivalent rate, and again below to turn
            // it into a per-timestep probability. The double conversion silently suppressed the
            // MM division probability far below its intended value. Now it's converted to a
            // step probability exactly once, matching how every other rate in the model
            // (e.g. T_CELL_DIV_RATE) is used.
            double pdiv = (1.0 / 1440) * MinToHour * scaleFactor;

            ///////////
            //MM dies//
            ///////////
            if (rn_BirthDeath < ProbScale(pdeath, TIMESTEP_AGENT)) {
                color = WHITE;
                Dispose();
            } else if (rn_BirthDeath < ProbScale(pdiv, TIMESTEP_AGENT)) {
                //////////////
                //MM divides//
                //////////////
                int emptyNeighbors = MapEmptyHood(MOORE_HOOD);
                if (emptyNeighbors > 0) {
                    int rIndex = G.rn.Int(emptyNeighbors);
                    int chosenCell = MOORE_HOOD[rIndex];
                    if (G.GetAgent(chosenCell) == null) {
                        simpleBoneCell child = G.NewAgentSQ(chosenCell);
                        child.type = this.type;
                        child.RESISTANT = this.RESISTANT;
                        child.TGFB_on = false;

                        // Assign this child a fresh unique clone ID, and record this cell as its
                        // parent. (Previously this incorrectly mutated the *parent's* cloneID in
                        // place and gave the child the same value, so IDs weren't stable/unique.)
                        child.parentID = this.cloneID;
                        child.cloneID = G.nextCloneID++;

                        // MHC-I loss
                        if (G.rn.Double() < G.antigenLoss && !this.mhcLoss) {
                            child.mhcIExpression = 0;
                            child.mhcLoss = true;
                        } else {
                            child.mhcIExpression = this.mhcIExpression;
                            child.mhcLoss = this.mhcLoss;
                        }
                        // BCMA loss
                        if (G.rn.Double() < G.antigenLoss && !this.bcmaLoss) {
                            child.bcmaExpression = 0;
                            child.bcmaLoss = true;
                        } else {
                            child.bcmaExpression = this.bcmaExpression;
                            child.bcmaLoss = this.bcmaLoss;
                        }
                        // GPRC5D loss
                        if (G.rn.Double() < G.antigenLoss && !this.gprc5dLoss) {
                            child.gprc5dExpression = 0;
                            child.gprc5dLoss = true;
                        } else {
                            child.gprc5dExpression = this.gprc5dExpression;
                            child.gprc5dLoss = this.gprc5dLoss;
                        }
                    }
                }
            }
        }

        if (type == naiveTcell) {
            this.myeloma_bound = false;
            this.lifespan = G.boundedGaussian(30, 1, 30, 34);

            if (time % (int) G.TIMESTEPS_PER_DAY == 0) {
                this.tcellAge += 1;
            }
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }

            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(MOORE_HOOD);
            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) &&
                    this.timeSinceDivision >= 12) {
                int chosenCell = MOORE_HOOD[G.rn.Int(emptyNeighbors)];
                if (G.GetAgent(chosenCell) == null) {
                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
                }
            }
            if (G.TCE) {
                if (G.rn.Double() > .5) {
                    this.pd_l1 = G.boundedGaussian(30, 1, 10, 40);
                    this.type = activeTcell;
                    return;
                }
            }
            for (int run = 0; run < 3; run++) {
                int moveToIndex = noisyseekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == extnaiveTcell) {
            this.myeloma_bound = false;
            this.lifespan = G.boundedGaussian(30, 1, 30, 34);

            if (time % (int) G.TIMESTEPS_PER_DAY == 0) {
                this.tcellAge += 1;
            }
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }

            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(MOORE_HOOD);
            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) &&
                    this.timeSinceDivision >= 12) {
                int chosenCell = MOORE_HOOD[G.rn.Int(emptyNeighbors)];
                if (G.GetAgent(chosenCell) == null) {
                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
                }
            }
            if (G.TCE) {
                if (G.rn.Double() > .5) {
                    this.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
                    this.type = EXHT_CELL;
                    return;
                }
            }
            for (int run = 0; run < 3; run++) {
                int moveToIndex = noisyseekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == activeTcell) {
            this.myeloma_bound = false;

            if (time % (int) G.TIMESTEPS_PER_DAY == 0) {
                this.tcellAge += 1;
            }
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }
            if (G.TGFB.Get(Isq()) >= (G.TGFBthresh * 2)) {
                this.pd_1 += 1;
            }
            // Exhaustion check
            if (this.pd_1 > pd_l1) {
                this.type = EXHT_CELL;
                return;
            }
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(MOORE_HOOD);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) &&
                    this.timeSinceDivision >= 6) {
                int chosenCell = MOORE_HOOD[G.rn.Int(emptyNeighbors)];
                if (G.GetAgent(chosenCell) == null) {
                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.pd_l1 = G.boundedGaussian(20, 1, 10, 40);
                }
            }

            for (int run = 0; run < 3; run++) {
                int[] movdivHood = MooreHood(true);
                int options = MapOccupiedHood(movdivHood);

                // Killing logic
                for (int j = 0; j < options; j++) {
                    if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM) {
                        if (G.BCMA_TCE && G.GetAgent(movdivHood[j]).bcmaExpression > 0
                                || G.GPRC5D_TCE && G.GetAgent(movdivHood[j]).gprc5dExpression > 0) {

                            double killProb = p_kill + (p_kill * G.GetAgent(movdivHood[j]).mhcIExpression);
                            if (G.rn.Double() < killProb) {
                                G.GetAgent(movdivHood[j]).Tcell_Kill();
                                this.pd_1 += 1;
                                this.PERF_on = true;
                            } else {
                                this.PERF_on = false;
                            }
                        }
                    }
                }

                int moveToIndex = seekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == EXHT_CELL) {
            double EXH_RECOG_NOISE = 0.5;

            // Restored: age-based death for exhausted T cells was previously dead code
            // (the increment was commented out), meaning an exhausted cell's `tcellAge` froze
            // at whatever value it had when it transitioned from activeTcell, and the
            // `tcellAge >= lifespan` death check below could never fire from natural aging.
            // If you want exhausted cells to *not* age out this way (e.g. only die via
            // pd_1 >= pd_l1 exhaustion), remove this block again.
            if (time % (int) G.TIMESTEPS_PER_DAY == 0) {
                this.tcellAge += 1;
            }
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }

            if (this.pd_1 >= this.pd_l1) {
                this.Dispose();
                return;
            }
            this.timeSinceDivision += TIMESTEP_AGENT;
            int emptyNeighbors = MapEmptyHood(MOORE_HOOD);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) &&
                    this.timeSinceDivision >= 24) {
                int chosenCell = MOORE_HOOD[G.rn.Int(emptyNeighbors)];
                if (G.GetAgent(chosenCell) == null) {
                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
                }
            }
            for (int run = 0; run < 3; run++) {
                int[] movdivHood = MooreHood(true);
                int options = MapOccupiedHood(movdivHood);
                for (int j = 0; j < options; j++) {
                    if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == MM) {
                        if (G.BCMA_TCE && G.GetAgent(movdivHood[j]).bcmaExpression > 0
                                || G.GPRC5D_TCE && G.GetAgent(movdivHood[j]).gprc5dExpression > 0) {
                            double killProb = extp_kill + (extp_kill * G.GetAgent(movdivHood[j]).mhcIExpression);
                            if (G.rn.Double() >= EXH_RECOG_NOISE) {
                                if (G.rn.Double() < killProb) {
                                    G.GetAgent(movdivHood[j]).Tcell_Kill();
                                    this.pd_1 += 1;
                                    this.PERF_on = true;
                                } else {
                                    this.PERF_on = false;
                                }
                            }
                        }
                    }
                }
                int moveToIndex = noisyseekCXCL9();
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }

        if (type == supressorTcell) {
            int[] hood = MooreHood(true);

            if (time % (int) G.TIMESTEPS_PER_DAY == 0) {
                this.tcellAge += 1;
            }
            if (this.tcellAge >= this.lifespan) {
                this.Dispose();
                return;
            }
            if (this.pd_1 >= this.pd_l1) {
                this.Dispose();
                return;
            }
            int emptyNeighbors = MapEmptyHood(MOORE_HOOD);

            if (emptyNeighbors > 0 &&
                    G.rn.Double() < ProbScale(G.T_CELL_DIV_RATE, TIMESTEP_AGENT) &&
                    this.timeSinceDivision >= 24) {

                int chosenCell = MOORE_HOOD[G.rn.Int(emptyNeighbors)];

                if (G.GetAgent(chosenCell) == null) {
                    simpleBoneCell child = G.NewAgentSQ(chosenCell);
                    child.type = this.type;
                    child.tcellAge = 0;
                    child.pd_l1 = G.boundedGaussian(10, 1, 10, 20);
                }
            }

            // Check neighbors for T cells
            int occupied = MapOccupiedHood(hood);
            simpleBoneCell target = null;
            for (int j = 0; j < occupied; j++) {
                simpleBoneCell neighbor = G.GetAgent(hood[j]);
                if (neighbor.type == activeTcell ||
                        neighbor.type == EXHT_CELL ||
                        neighbor.type == extnaiveTcell ||
                        neighbor.type == naiveTcell) {
                    target = neighbor;
                    break;
                }
            }

            TGFB_on = target != null;

            int moveTo = seekPerf();
            if (moveTo >= 0 && G.GetAgent(moveTo) == null) {
                MoveSQ(moveTo);
            }
        }
    }
}
