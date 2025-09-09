package boneRemodeling_2022May17;
import HAL.GridsAndAgents.*;
import HAL.Gui.*;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;
import static boneRemodeling_2022May17.BoneGrid_2022May17.*;
import static HAL.Util.*;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     GRID CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public class BoneGrid_2022May17 extends AgentGrid2D<BoneCell_2022May17> implements SerializableModel {

    ///////////////
    //GRID FIELDS//
    ///////////////

    //3/4/23 testing pob differentiation at 0.004 for 25 sims
    //4.19/24 need to check the increase in aOC resorb. Now it is 5, running without myeloma and setting it to 0.1
    //4/23/24 testing new aOC aging effect, increased aOC lifespan - this is done 4/24/24
    //4/24/24 running baes case myeloma seeded at 50 years
    // myeloma seeded at year 10 still on the cluster
    // 4/25/24 runing base case with OC increased lifespan, no myeloma
    // running MSC protection on MM with osteolcast lifespan increased
    // 4/28/24 reran baseline aged 85 years, no myeloma, decrease pOB proliferation, increase aOC lifespan
    // 4/28/24 rerunning msc benefit, aged 85 years, no myeloma, decrease pOB proliferation, increase aOC lifespan
    // 4/29/24 running myeloma seeded at year 30 in the model, no msc benefit, no oc benefit, only decreased pOB proliferation
    // 4/29/24 added T-cells
    // 4/30/24 Need to run myeloma with msc benefit at model timepoint 50, no OC benefit, decreased pOB proliferation, running in myeloma 3
    // 4/30/24 need to run same as above but with no msc benefit, running in myeloma2 folder
    // 5/1/24 added = to the myeloma seed time value to test tcells
    // 5/7/24 testing Tcells for 25 iterations for 1 year, running baseline for 1 year, both with myeloma
    // baseline with myeloma is in myeloma2
    // 5/7/24 need to run the condition where there is an msc benefit to myeloma over again, in myeloma3
    // 5/8/24 rerunning base condition with msc benefit to myeloma over again in myeloma3, seed time was wrong on the last run
    // 5/8/24 rerunning Tcells with different pd-1 value
    // 5/8/24 Need to rerun the 1 year myeloma over again, unintended 0 value in msc benefit caused year 0 to have extrmeely low level of proliferation
    // 5/13/24 changed to Tcells to move 30um per timestep, running experiments now
    // 5/13/24 added two null checks, use comand f to find them if they need to be changed later
    // 5/17/24 running aged myeloma with Tcells for 85 years, more so to test if the cluster is working
    // 5/17/24 running the above test in myeloma2 cluster folder
    // 5/17/24 modified Tcell code to exit outer loo if a myeloma cell is encountered
    // 5/20/24 running myeloma for 4 years with and without Tcells in Tcells and myeloma2 respectively
    // 10/24/24 ran antigen loss rate of 10^-3
    //10/24/24 ran antigen loss rate of 10^-5
    //10/24/24 ran antigen loss rate of 10^-7
    //10/24/24 need to rerun the baseline condition with no Tcells and antigen loss rate of 0
    //1/22/25 added a feature to prevent infinite loops in CD8 T-cell placement
    //1/23/25 Changed the number of cancer cells CD8 T-cells can kill before exhuastion to ~100



    //Variables to switch on/off treatment/other things
    public static boolean TGFB_INHIBITOR = false;
    public static boolean RANKL_INHIBITOR = false;
    public static boolean BISPHOSPHONATE = false;
    public static boolean BORTEZOMIB = false;
    public static boolean MYELOMA = true;
    public static boolean TCELL = true;
    public static boolean EMDR = false;
    public static boolean TREATMENT_ON = false; //this is to control treatment on/off timer in MAIN

    //CLUSTER
    public static boolean PARAM_SWEEP = false; //use when importing parameters to loop through
    public static boolean HEADLESS = false; //use true with cluster
    public static boolean LOCAL = true; // use false with cluster
    public static double numSteps = 2.0*365.0*24.0*60.0; // years the model will run
    public static int numSims = 25; //Number of Simulations
    public final static int BONE = RGB256(255,255,250), MSC = RGB256(135,206,250),
            pOB = RGB256(100,149,237), aOB = BLUE, pOC = RGB256(230,100,130),
            aOC = RED, LINING = RGB256(64,106,151), MM = RGB256(0,128,0),
            activeTcell = RGB256(17, 150, 150),
            EXHT_CELL=RGB256(200, 50, 250),
            supressorTcell =RGB256(255, 165, 0),
            bloodVessel=RGB256(138, 3, 3),
            ALMM=RGB(255, 105, 180);


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
    public PDEGrid2D IFNG;
    public int MM_Division_List_Size = 10;//10 timesteps = 60 min = 1 hr;

    public int[] tmoveHood = MooreHood(true);
    public ArrayList<Integer> InitBoneList = new ArrayList<>();
    public ArrayList<BoneCell_2022May17> AllBoneList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<BoneCell_2022May17> LiningList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<BoneCell_2022May17> tempEventList = new ArrayList<>(); //This list temporarily stores the nearest bone-lining cell neighbors

    boolean [][] exposedBone = new boolean[xDim][yDim];

    //public int eOpt;
    FileIO out;
    FileIO InitialBone;
//    FileIO params;

    //This is important for serializable model
    @Override
    public void SetupConstructors(){
        this._PassAgentConstructor(BoneCell_2022May17.class);
    }

    ////////////////////
    //GRID CONSTRUCTOR//
    ////////////////////

    public BoneGrid_2022May17(int xDim, int yDim, Rand rn, String Bone_FileName) {
        super(xDim, yDim, BoneCell_2022May17.class,true,true);
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
    //5. CollectLINING     //
    //6. Draw              //
    //7. DrawRANKL         //
    //8. RecordRANKL       //
    //9. RecordOut         //
    /////////////////////////


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


        if(mode=="w") {

            out.Write("Timestep" + "," + "BONE" + "," + "pOB" + "," + "aOB" + "," + "pOC" + "," + "aOC" + "," + "MSC" + "," + "LINING" + "," + "S_MM" + "," + "R_MM"+"," +"AL_MM"+ ","+"TCell"+","+"ExtTcell"+"," +"T-reg" +","+ "TREATMENT_ON" + "," + "BORTEZOMIB" + "," + "MYELOMA" + "\n");
        }

    }

    public void closeFileIO () {

        out.Close();
//        params.Close();
    }

    public void SetParams(int prow, ArrayList<String> param_list){
        //returns an array list of all lines from the file as stringsftype == lining

        String[] split_param_list = param_list.get(prow).split(",");


        MYELOMA = Boolean.parseBoolean(split_param_list[0]);
        BORTEZOMIB = Boolean.parseBoolean(split_param_list[1]);
        pmutate = Double.parseDouble(split_param_list[2]);
        EMDR = Boolean.parseBoolean(split_param_list[3]);
        dose = Double.parseDouble(split_param_list[4]);
    }

    // Helper method to check if a cell is in the bone lining area
    private boolean isInLining(int x, int y) {
        BoneCell_2022May17 cell = GetAgent(x, y);
        return (cell != null && cell.type == LINING);
    }
    public void InitBone() {


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

        //int myelomaCellsToPlace = 3000; // Total number of myeloma cells to place
        int myelomaCellsToPlace = 3000; // Total number of myeloma cells to place
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
                        GetAgent(xNeighbor, yNeighbor).mhc_i_expression = boundedGaussian(0.5, 0.1,0, 1);
                        if (rn.Double() < bcmaNegfraction){
                            GetAgent(xNeighbor, yNeighbor).bcmaLoss = true;

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

    public void ModelStep(int time, double [] Cell_Counts) {

        //STEP 0: UPDATE GRIDTICK
        /////////////////////////////////////////////////
        //STEP 1: REACTION-DIFFUSION EQUATION FOR RANKL//
        /////////////////////////////////////////////////
        int i=0;
        double stol = 1.0e-6;//1e-6; //steady-state tolerance


        for (int x = 0; x < CXCL9.xDim; x++) {
            for (int y = 0; y < CXCL9.yDim; y++) {
                if (GetAgent(x,y)!=null && GetAgent(x, y).type == MM && GetAgent(x, y).bcmaLoss!=true )  {
                    CXCL9.Add(x, y, CXCL9_productionRate/maxCXCL9 );
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
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
                    if (GetAgent(x, y).myeloma_bound==true){
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
        for (BoneCell_2022May17 c: this) {
            c.CellStep(time, Cell_Counts);
        }

    }

    ////////////////////////////
    //Full Domain (No zoom-in)//
    ////////////////////////////
    public void Draw(UIGrid vis, UILabel days, int i) {
        days.SetText("days: "+i*convert_to_days);
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                if (drawMe != null && drawMe.type==MM && drawMe.RESISTANT==true){
                    vis.SetPix(x,y, BLACK);
                } else if (drawMe != null && drawMe.type==MM && drawMe.bcmaLoss==true){
                    vis.SetPix(x,y, BLACK);
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
                BoneCell_2022May17 drawMe = GetAgent(x,y);
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
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else
                    vis.SetPix(x, y, HeatMapRGB(IFNG.Get(x, y)));

            }
        }
    }


    public void RecordOut(FileIO writeHere,int time, boolean treatment_on, boolean btz, boolean myeloma){
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int[] cts = new int[13];

        for (BoneCell_2022May17 c : this) {
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

        }
        //population of one timestep per line
        //writeHere.Write(ct_BONE+","+ct_pOB+","+ct_aOB+","+ct_pOC+","+ct_aOC+","+ct_MSC+","+ct_LINING+"\n");
        writeHere.Write(time+",");
        writeHere.WriteDelimit(cts,",");
        writeHere.Write("," + treatment_on + "," + btz + "," + myeloma + "\n");
    }

    public double[] CellCounts(){
        //int BONE = 0, pOB = 1, aOB = 2, pOC = 3, aOC = 4, MSC = 5, LINING = 6, MM = 7;
        double[] cts = new double[13];

        for (BoneCell_2022May17 c : this) {
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

    public static void main(String[] args) {

//        boolean HEADLESS = true; //use true with cluster

        if (HEADLESS) {
            System.setProperty("java.awt.headless", "true");
        }

        int xDim = 160;//500; //px
        int yDim = 150;//250; //px
        String sdf = new SimpleDateFormat("yyyyMMMd").format(new Date());

        UIWindow win = HEADLESS ? null : new UIWindow("Normal Bone Remodeling");
        String fn = "Bone_" + sdf;
        File dir = new File(fn);
        dir.mkdir();

        ///////////////////////
        //FOR PARAMETER SWEEP//
        ///////////////////////

        int param_list_size;
        ArrayList<String> param_list = null;

        if(PARAM_SWEEP) {
            FileIO Params = new FileIO("Bone/boneRemodeling_2022May17/params.csv", "r");
            param_list = Params.Read();
            param_list_size = param_list.size();
        }else{
            param_list_size = 2; //this will go through iteration once for pre-defined parameters
        }


        //Define parameters
        for (int prow=1; prow<param_list_size; prow++) {

            for (int sim = 0; sim < numSims; sim++) {
                String subfolder = fn + "/Sim" + sim + "_row" + prow + "/";
                dir = new File(subfolder);
                boolean success = dir.mkdir();

                if (success) {


                    ////////////////////////////////////////
                    //UIGrid                              //
                    //compx=number of columns (default 1) //
                    //compy=number of rows (default 1)    //
                    ////////////////////////////////////////

                    //Full domain
                    UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
                    UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
                    UIGrid IFNG_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2


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

                    // GIF MAKER
                    String projPath = subfolder; //PWD() + subfolder;
                    GifMaker gm_Cell_vis = new GifMaker(projPath.concat("/").concat("CellVid").concat(".gif"), 100, true);
                    GifMaker gm_CXCL9_vis = new GifMaker(projPath.concat("/").concat("CXCL9").concat(".gif"), 100, true);
                    GifMaker gm_IFNG_vis = new GifMaker(projPath.concat("/").concat("IFNG").concat(".gif"), 100, true);
                    String Bone_Filename = null;
                    if (LOCAL) {
                        Bone_Filename = "/Users/80024703/Desktop/code/Bone/BAout_2020May5_Sim14.csv";
                    } else {
                        Bone_Filename = "Bone/BAout_2020May5_Sim14.csv";
                    }
                    BoneGrid_2022May17 g = new BoneGrid_2022May17(xDim, yDim, new Rand(), Bone_Filename); //set seed to reproduce results

                    //Rand(seed:1) when want to reproduce results

                    //Set treatment schedule
                    ArrayList<Double> Tx_subStart = new ArrayList<>();
                    if (g.Tx_Duration == g.Tx_Interval) { //continuous treatment
                        Tx_subStart.add(0.0);//1440.0/TIMESTEP_AGENT); //1 day=1440 minutes
                    } else { //pulsed treatment
                        Tx_subStart.add(0.0);//1440.0/TIMESTEP_AGENT); //1 day=1440 minutes
                    }
                    double subStart_time = 0;


                    //Set parameters
                    if (PARAM_SWEEP) { //Note: Other parameters dependent on this list will not be updated
                        g.SetParams(prow, param_list);
                    }


                    //Record Output
                    g.newFileIO(projPath, "w");

                    //Initialize model
                    g.InitBone();


//                  RMevents=new int[Nevents];
//                  g.rn.RandomIS(RMevents, 0, g.curI * g.TURNOVER_TIME);

                    List<Integer> vessLocations = new ArrayList<>();


                    double reducedMeanFraction = 0.23; // 0.23 (23%)is the default  based on experimental data
                    double tregFraction = .11; //.18 is the default based on experimental data 4/7/25 was set to .11. 4/17/25 .11 is now the default
                    double tcellPercent = 0.868/100 * (xDim * yDim);
                    double tregPercent = (((0.6 * tregFraction ) / 100)  * (xDim * yDim));
                    int activeTcells = (int) Math.floor((1-reducedMeanFraction) * tcellPercent);
                    int activeExtcells = (int) Math.floor(reducedMeanFraction * tcellPercent); //208 total Tcells in the model
                    int activeTregs = (int) Math.floor(tregPercent);
                    double dailyExhaustionProb = 0.00131;
                    double dailyTcellDecreaseProb = 0.000583;
                    double dailyTregIncreasenProb = 0.000791;
                    int tcellDecrease = 208;
                    //Loop through//
                    ////////////////
                    boolean initial_recruitment = false;
                    boolean initial_vessel = false;
                    for (int i = 0; i < g.Nts; i++) {//g.Nts
                        timeStep = i; // 240 ts is equal to one day

                        if (!HEADLESS) {
                            win.TickPause(10); //slow down visualization
                        }

                        double[] Cell_Counts = g.CellCounts(); //int BONE = 0, pOB = 1, aOB = 2, pOC = 3, aOC = 4, MSC = 5, LINING = 6, MM = 7;
                        if (initial_vessel==false) {
                            int vesselNumber = 200;
                            int k = 0;
                            while (k < vesselNumber) {
                                int xinit = g.rn.Int(xDim);
                                int yinit = g.rn.Int(yDim);
                                if (g.GetAgent(xinit, yinit) != null && g.GetAgent(xinit, yinit).type != LINING&& g.GetAgent(xinit, yinit).type != BONE) {
                                    g.GetAgent(xinit, yinit).Dispose();
                                }
                                while (g.PopAt(xinit, yinit) > 0) {
                                    xinit = g.rn.Int(xDim);
                                    yinit = g.rn.Int(yDim);
                                }
                                // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                                BoneCell_2022May17 c = g.NewAgentSQ(xinit, yinit);
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
                                            BoneCell_2022May17 c = g.NewAgentSQ(chosenCell);
                                            c.type = activeTcell;
                                            if (g.rn.Double() < reducedMeanFraction) {
                                                c.type = EXHT_CELL;
                                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            } else {
                                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
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
                                            BoneCell_2022May17 c = g.NewAgentSQ(candidate);
                                            c.type = activeTcell;
                                            if (g.rn.Double() < reducedMeanFraction) {
                                                c.type = EXHT_CELL;
                                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            } else {
                                                c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            }
                                            j++;
                                            break;
                                        }
                                        tries++;
                                    }
                                }

                                if (attempts >= maxAttempts) {
                                    System.out.println("Reached maximum attempts. T-cells placed: " + j + " out of " + InitTcells);
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
                                            BoneCell_2022May17 c = g.NewAgentSQ(chosenCell);
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
                                            BoneCell_2022May17 c = g.NewAgentSQ(candidate);
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
                            if (Cell_Counts[10] < activeTcells) {
                                InitTcells = (int) Math.floor(g.boundedGaussian((activeTcells - Cell_Counts[10]), 4, activeTcells - Cell_Counts[10]-4, activeTcells - Cell_Counts[10] + 1));
                            } else {
                                //InitTcells = (int) g.boundedGaussian(1, 1, 1, 10);
                                InitTcells = 0;
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

                                            BoneCell_2022May17 c = g.NewAgentSQ(chosenCell);
                                            c.type = activeTcell;
                                            c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            j++;
                                            placed = true;
                                            num_placed+= 1;
                                        }
                                    }
                                }

                                // Fallback to boundary if not placed
                                if (!placed) {
                                    int tries = 0;
                                    while (tries < boundaries.length) {
                                        int randomBoundaryIdx = rand.nextInt(boundaries.length);
                                        int candidate = boundaries[randomBoundaryIdx];
                                        if (g.GetAgent(candidate) == null) {
                                            BoneCell_2022May17 c = g.NewAgentSQ(candidate);
                                            c.type = activeTcell;
                                            c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            num_placed+=1;
                                            j++;
                                            break;
                                        }
                                        tries++;
                                    }
                                }
                            }
//                            System.out.println("Active Tcells");
//                            System.out.println(activeTcells);
//                            System.out.println(Cell_Counts[10]);
//                            System.out.println("Placed Tcells "+num_placed);
//                            System.out.println();

                        }

                        if (TCELL && initial_recruitment && daysPassed > 0) {
                            int InitTcells;
                            if (Cell_Counts[11] < activeExtcells) {
                                InitTcells = (int) Math.floor(g.boundedGaussian((activeExtcells - Cell_Counts[11]), 3, activeExtcells - Cell_Counts[11]-3, activeExtcells - Cell_Counts[11] + 1));
                            } else {
                                //InitTcells = (int) g.boundedGaussian(1, 1, 1, 10);
                                InitTcells = 0;
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

                                // Try blood vessel neighborhood
                                int randomVessLocation = vessLocations.get(rand.nextInt(vessLocations.size()));

                                if (g.GetAgent(randomVessLocation).type == bloodVessel) {
                                    int[] movdivHood = MooreHood(true);
                                    int emptyNeighbors = g.MapEmptyHood(movdivHood, randomVessLocation);

                                    if (emptyNeighbors > 0) {
                                        int chosenIndex = g.rn.Int(emptyNeighbors);
                                        int chosenCell = movdivHood[chosenIndex];

                                        if (g.GetAgent(chosenCell) == null) {
                                            BoneCell_2022May17 c = g.NewAgentSQ(chosenCell);
                                            c.type = EXHT_CELL;
                                            c.pd_1 = g.boundedGaussian(10, 1, 1, 20); // Lower mean for exhausted
                                            num_placed+=1;
                                            j++;
                                            placed = true;
                                        }
                                    }
                                }

                                // Fallback to random empty boundary site if no vessel placement succeeded
                                if (!placed) {
                                    int tries = 0;
                                    while (tries < boundaries.length) {
                                        int randomBoundaryIdx = rand.nextInt(boundaries.length);
                                        int candidate = boundaries[randomBoundaryIdx];
                                        if (g.GetAgent(candidate) == null) {
                                            BoneCell_2022May17 c = g.NewAgentSQ(candidate);
                                            c.type = EXHT_CELL;
                                            c.pd_1 = g.boundedGaussian(10, 1, 1, 20);
                                            num_placed+=1;
                                            j++;
                                            break;
                                        }
                                        tries++;
                                    }
                                }

                            }
//                            System.out.println("Ext Tcells");
//                            System.out.println(activeExtcells);
//                            System.out.println(Cell_Counts[11]);
//                            System.out.println("Placed Ext Tcells "+num_placed);
//                            System.out.println();
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
                                            BoneCell_2022May17 c = g.NewAgentSQ(chosenCell);
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
                                        int candidate = boundaries[randomBoundaryIdx];
                                        if (g.GetAgent(candidate) == null) {
                                            BoneCell_2022May17 c = g.NewAgentSQ(candidate);
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
//                            System.out.println("Suppresor Tcells");
//                            System.out.println(activeTregs);
//                            System.out.println(Cell_Counts[12]);
//                            System.out.println("Placed Supp Tcells "+ num_placed);
//                            System.out.println();
                        }

                        if (g.rn.Double() <dailyTcellDecreaseProb ){
                            tcellDecrease = tcellDecrease -1;
                            activeTcells = activeTcells-1;
                            activeExtcells = activeExtcells -1;
                        }
                        if ( g.rn.Double() < dailyTregIncreasenProb) {
                            activeTregs = activeTregs + 1;
                        }
                        if (activeExtcells <= tcellDecrease && g.rn.Double() < dailyExhaustionProb) {
                            activeTcells = activeTcells - 1;
                            activeExtcells = activeExtcells + 1;
                        }

                        g.ModelStep(i, Cell_Counts);

                        g.Draw(Cell_vis, days, i);
                        g.DrawCXCL9(CXCL9_vis);
                        g.DrawIFNG(IFNG_vis);

                        if (i % (24.0 * 60.0 / (MinToHour * TIMESTEP_AGENT)) == 0) {
                            //240ts = 1 day



                            g.RecordOut(g.out, i, TREATMENT_ON, BORTEZOMIB, MYELOMA);
                            //System.out.println("Timestep =" + i);
                            daysPassed += 1;
                        }
                    }

                    //Close FileIO
                    g.closeFileIO();

                    //Close GifMaker
                    gm_Cell_vis.Close();
                    gm_CXCL9_vis.Close();

                    //Close UIWindow
                    if (!HEADLESS) {
                        win.Close();
                    }
                }
            } //END OF MULTIPLE RUNS
        } //end parameter loop
    }


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     CELL CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BoneCell_2022May17 extends AgentSQ2Dunstackable<BoneGrid_2022May17> {

    ///////////////
    //CELL FIELDS//
    ///////////////

}
