package cd8Test;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.*;
import HAL.Rand;
import HAL.Tools.FileIO;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import static cd8Test.gridLattice.*;
import static HAL.Util.*;

//////////////
//GRID CLASS//
//////////////

public class gridLattice extends AgentGrid2D<agentLattice> {
    public static int timeStep = 0; // initial timestep
    public final static double STEPS = (28 * 24.0); // timesteps in hours
    public final static int ITERATIONS =1 ; // number of times the model will run
    public static boolean VISUALS = true; // turns on/off visualizations for the model
    public final static int activeTcell = 1; // active cd8 t-cells
    public final static int Myeloma = 2; // myeloma cell
    public final static int BONE = 3; // bone
    public final static int exhaustedTcells = 4; // exhausted cd8 t-cell
    public final static int LINING = 5; // bone lining
    public final static int inactiveTcell = 6; // inactive tcells
    public final int InitactiveTcell = 15; // initial number of active tcells
    public final int InitinactiveTcells = 0; //initial number of inactive tcells
    public final int InitexhaustedTcells = 5; // initial number of exhausted tcells
    final static double Myeloma_PROLIFERATION_PROB = 1.0 /24; //cancer cell division rate
    final static double Myeloma_DEATH_PROB = 1.0 /72; // cancer cell apoptosis rate
    final static double mutProb = 0.0000007; // probability that a myeloma cell will lose BCMA due to a mutation
    public double T_CELL_DEATH_RATE = 1.0 /24; // tcell death rate - used for exhausted tcells
    double CXCL9_productionRate = .1; // production of CXCL9
    double CXCL9_decayRate = -0.3; // decay of CXCL9
    double CXCL9_DiffCoef = 2.5; // diffusion coefficient for CXCL9
    double maxCXCL9 = 1.7; // used for numerical stability
    double Tcell_TaxisCoeff = 3.0e9; // tcell movement towards CXCL9
    double Tcell_DiffCoef = 0.1*3.0; // diffusion coefficient for cd8 t-cells

    boolean TCE_ON = false; // switching on tce therapy

    boolean preRes = false;

    double resFrac = 0;

    double maxTCE = 0; // not using this yet

    double TCE_start = 0; // not using this yet

    double TCE_Duration = 0; // not using this yet

    public PDEGrid2D CXCL9; // pde grid for CXCL9

    public PDEGrid2D TCE; // not using this yet

    public int[] tmoveHood = MooreHood(true); //Have option of no movement


    public Rand rn = new Rand();

    FileIO output;

    ////////////////////
    //GRID CONSTRUCTOR//
    ////////////////////
    public double boundedGaussian(double mean, double dev, double min, double max) {
        double gauss = rn.Gaussian(0, 1);
        double val = dev * gauss + mean;
        while (val > max || val < min) {
            gauss = rn.Gaussian(0, 1);
            val = dev * gauss + mean;
        }
        return val;
    }

    boolean isWithinDistanceFromBone(int x, int y, int startX, int startY, int sideLength, int distance) {
        int boneEndX = startX + sideLength;
        int boneEndY = startY + sideLength;

        int minX = Math.max(0, startX - distance);
        int maxX = Math.min(xDim - 1, boneEndX + distance);
        int minY = Math.max(0, startY - distance);
        int maxY = Math.min(yDim - 1, boneEndY + distance);

        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }

    public gridLattice(int xDim, int yDim, String outFileName) {
        super(xDim, yDim, agentLattice.class, true, true);
        CXCL9 = new PDEGrid2D(xDim, yDim,true,true); //This assumes PERIODIC BOUNDARY CONDITION (wrapX=TRUE,wrapY=TRUE)
        output = new FileIO(outFileName, "w");
        output.Write("Time,TCELLS, EXTCELLS, MYELOMA, Inactive Tcells, CXCL9" + "\n");

        double totalArea = xDim * yDim;
        double blackBoxArea = 0.129 * totalArea;
        int sideLength = (int)Math.sqrt(blackBoxArea);

        // Calculate the center of the grid
        int centerX = xDim / 2;
        int centerY = yDim / 2;

        // Calculate the top-left corner of the black box
        int startX = centerX - sideLength / 2;
        int startY = centerY - sideLength / 2;

        // Place the black box border and inside color
        for (int i = startX; i < startX + sideLength; i++) {
            for (int j = startY; j < startY + sideLength; j++) {
                if (i == startX || i == startX + sideLength - 1 || j == startY || j == startY + sideLength - 1) {
                    // Set the cell type to BLACK_BOX and color to black
                    agentLattice c = NewAgentSQ(i, j);
                    c.type = LINING;
                } else {
                    agentLattice c = NewAgentSQ(i, j);
                    c.type = BONE;
                }
            }
        }
        for (int i = 0; i < 500; i++) {
            // Recruit one myeloma cell within 10 grid spaces from the bone
            int myelomaX, myelomaY;
            do {
                myelomaX = rn.Int(xDim);
                myelomaY = rn.Int(yDim);
            } while (!isWithinDistanceFromBone(myelomaX, myelomaY, startX, startY, sideLength, 5) || PopAt(myelomaX, myelomaY) > 0);

            agentLattice c = NewAgentSQ(myelomaX, myelomaY);
            c.type = Myeloma;
            c.SetCellColor();
        }
    }

    ////////////////
    //GRID METHODS//
    ////////////////
    public void Draw(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                agentLattice drawMe = GetAgent(x, y);
                if (drawMe != null) {
                    vis.SetPix(x, y, drawMe.color);
                } else {
                    vis.SetPix(x,y, RGB256(240, 220, 220));
                }
            }
        }
    }
    public void DrawCXCL9(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                agentLattice drawMe = GetAgent(x, y);
                if (drawMe != null) {
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));
                } else {
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));

                }
            }
        }
    }

    public double[] CellCounts() {
        double[] counts = new double[5];
        for (agentLattice c : this) {
            if (c.type == activeTcell) {
                counts[0]++;
            }
            if (c.type == exhaustedTcells) {
                counts[1]++;
            }
            if (c.type == Myeloma) {
                counts[2]++;
            }
            if (c.type == inactiveTcell) {
                counts[3]++;
            }
        }
        return counts;
    }

    public void RecordOut(FileIO writeHere, int time, double[] cts) {
        writeHere.Write(time + ",");
        writeHere.WriteDelimit(cts, ",");
        writeHere.Write("," + "\n");
    }


    public void ModelStep(double timeSteep, double day) {
        // if ( day > TCE_start){
        //         if (timeStep % 24 == 0) {
        for (int x = 0; x < CXCL9.xDim; x++) {
            for (int y = 0; y < CXCL9.yDim; y++) {
                if (GetAgent(x,y)!=null && GetAgent(x, y).type == Myeloma) {
                    CXCL9.Add(x, y, CXCL9_productionRate/maxCXCL9 );

                }
            }
        }

        //TCE Diffusion
        CXCL9.DiffusionADI(CXCL9_DiffCoef);

        //Natural Decay of RANKL
        CXCL9.MulAll(CXCL9_decayRate);
        CXCL9.Update();

        ShuffleAgents(rn);
        for (agentLattice c : this) {
            c.CellStep(day);
        }
    }

    ////////
    //MAIN//
    ////////

    public static void main (String[]args) throws IOException {
        int j = 0;
        //TIME STAMP
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH_SS");
        java.util.Date now = new java.util.Date();
        String date_time = dateFormat.format(now);
        int numIterations = ITERATIONS; // Change this value to specify the number of iterations
        while (j < numIterations){
            double day=0;
            String save = new String();
            //GridWindow win = null;
            int xDim = 160;
            int yDim = 150;
            double TIMESTEPS = STEPS;

            // PATHS
            String projPath = "/Users/80024703/Desktop/tcellTest/";
            String setting_dir = projPath + save + "/output/" + "/" + date_time + "/";
            // CREATE OUTPUT DIR
            new File(setting_dir).mkdirs();
            String output_dir = setting_dir;
            new File(output_dir).mkdirs();
            // OUTPUT FILES
            String path_to_output_file = output_dir.concat("/").concat("CellCounts_" + j).concat(".csv");
            // OUTPUT WINDOW
            UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
            UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 4,2,5);
            UIWindow win = new UIWindow("");
            win.AddCol(0, new UILabel("Cells"));
            win.AddCol(0, Cell_vis);
            win.AddCol(2, new UILabel("CXCl9"));
            win.AddCol(2, CXCL9_vis);
            win.RunGui();
            GifMaker gm_Cell_vis = new GifMaker(output_dir.concat("/").concat("CellVid").concat(".gif"), 100, true);
            GifMaker gm_CXCL9_vis = new GifMaker(output_dir.concat("/").concat("CXCL9Vid").concat(".gif"), 100, true);
            // GRID
            gridLattice g = new gridLattice(xDim, yDim, path_to_output_file);

            boolean initial_recruitment = true;
            // TIME LOOP
            for (int i = 0; i < TIMESTEPS; i++) {;
                timeStep = i;
                if (VISUALS) {
                    win.TickPause(100); //slows or speeds up the video frame rate
                }
                // COUNT CURRENT POPULATION
                double[] cell_count = g.CellCounts();
                if (cell_count[2] >= 500) { //recruiting tcells
                    if (initial_recruitment) {
                        int InitTcells = g.InitactiveTcell;
                        int k = 0;
                        while (k < InitTcells) {
                            int xinit = g.rn.Int(xDim);
                            int yinit = g.rn.Int(yDim);
                            while (g.PopAt(xinit, yinit) > 0) {
                                xinit = g.rn.Int(xDim);
                                yinit = g.rn.Int(yDim);
                            }
                            // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                            agentLattice c = g.NewAgentSQ(xinit, yinit);
                            c.type = activeTcell;
                            c.SetCellColor();
                            k++;
                        }
                    }
                    initial_recruitment = false;
                }
                if(timeStep == 0){
                    int k = 0;
                    while (k < g.InitinactiveTcells) {
                        int xinit = g.rn.Int(xDim);
                        int yinit = g.rn.Int(yDim);
                        while (g.PopAt(xinit, yinit) > 0) {
                            xinit = g.rn.Int(xDim);
                            yinit = g.rn.Int(yDim);
                        }
                        // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                        agentLattice c = g.NewAgentSQ(xinit, yinit);
                        c.type = inactiveTcell;
                        c.SetCellColor();
                        k++;
                    }

                    int m = 0;
                    while (m < g.InitexhaustedTcells) {
                        int xinit = g.rn.Int(xDim);
                        int yinit = g.rn.Int(yDim);
                        while (g.PopAt(xinit, yinit) > 0) {
                            xinit = g.rn.Int(xDim);
                            yinit = g.rn.Int(yDim);
                        }
                        // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                        agentLattice c = g.NewAgentSQ(xinit, yinit);
                        c.type = exhaustedTcells;
                        c.SetCellColor();
                        m++;
                    }

                }
                cell_count[4] = g.CXCL9.GetAvg();
                // RUN MODEL SIMULATION STEP
                g.ModelStep(timeStep, day);
                // GRAPHICAL OUTPUT
                if (VISUALS) {
                    g.Draw(Cell_vis);
                    g.DrawCXCL9(CXCL9_vis);
                }
                // OTHER OUTPUT
                if (i % 24.0 == 0) {
                    day+= 1;
                    System.out.println("Day " + day);
                    g.RecordOut(g.output,i/24, cell_count);
                    //System.out.println("Day "+i/240);
                }
            }
            g.output.Close();
            gm_Cell_vis.Close();
            gm_CXCL9_vis.Close();
            win.Close();

            if (j < numIterations){
                System.out.println("Current Iteration: " + j);}
            j++;
        }
    }
}

//////////////
//CELL CLASS//
//////////////

class agentLattice extends AgentSQ2Dunstackable<gridLattice> {
    public int type; // type of cell
    int color; // color of each cell for visualization
    double pd_1 = 0; // pd_1 for tcells
    double pd_l1 = 0; // pd_l1 for tcells
    double tcellAge = 0; // used for tracking the age of a tcell
    double hours_since_kill =0; // hours since a tcell last killed a myeloma cell
    double hours_since_division = 0; //days since the Tcell has divided
    boolean BCMA = true; // used for whether or not a myeloma cell is expressing bcma

    ////////////////
    //CELL METHODS//
    ////////////////
    final public static int tcell = RGB256(17, 150, 150); //tcell color
    final public static int exhausted = RGB256(200, 50, 250); //exhausted tcell color
    final public static int cancer = RGB256(47, 32, 66); // myeloma cell color

    final public static int bone = RGB256(64,106,151); //bone color
    final public static int inactivetcell = RGB256(255, 165, 0);
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
    public static int CategorialColor(int index) {
        switch (index) {
            case activeTcell:
                return tcell;
            case exhaustedTcells:
                return exhausted;
            case Myeloma:
                return cancer;
            case BONE:
                return bone;
            case inactiveTcell:
                return inactivetcell;
            default:
                throw new IllegalArgumentException("index outside color category range index: " + index);
        }
    }
    void SetCellColor() {
        color = CategorialColor((int) Math.round(0 + type));
    }
    void Tcell_Kill() {
        if (type == Myeloma) {
            this.Dispose();}
    }
    void CellStep(double day) {

        if (type == BONE){
            this.color = RGB256(255,255,250);
        }

        if (type == LINING){
            this.color = RGB256(64,106,151);
        }

        if (type == exhaustedTcells) {
            this.SetCellColor();
            int[] movdivHood = MooreHood(true); // for division and movement
            int emptyNeighbors = MapEmptyHood(movdivHood); // mapping empty spots
            double rn_BirthDeath = G.rn.Double();
            double pdiv = G.T_CELL_DEATH_RATE;
            if (emptyNeighbors > 0) {
                int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                    MoveSQ(chosenCell);
                }
            }
            if (rn_BirthDeath < ProbScale(pdiv, 1/24)) {
                this.Dispose();
            }
        }

        if (type == activeTcell) {
            boolean encounteredMyeloma = false; // Track if a myeloma cell is encountered
            boolean agedDeath = false; // Track if the T cell is aged and should die

            // Initialize pd_1 value if the T cell is newly created
            if (this.tcellAge == 0) {
                this.pd_1 = G.boundedGaussian(10, 5, 5, 20);
            }

            // Increment T cell age and hours since last kill every hour
            this.tcellAge += 1;
            this.hours_since_kill += 1;
            this.hours_since_division += 1;

            for (int run = 0; run < 3; run++) {
                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots
                int emptyNeighbors = MapEmptyHood(movdivHood); // Mapping empty spots

                // Age-related death logic
                if (this.tcellAge == 720 && emptyNeighbors > 0) { // Assuming 720 hours is the lifespan (30 days)
                    for (int i = 0; i < emptyNeighbors; i++) {
                        int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                        int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                        if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                            agentLattice child = G.NewAgentSQ(chosenCell);
                            child.type = this.type;
                            child.SetCellColor();
                            child.tcellAge = 0;
                            child.pd_1 = Math.floor(this.pd_1 / 4);
                            child.pd_l1 = this.pd_l1;
                            child.hours_since_division = 0;
                            break;
                        }
                    }
                    this.Dispose();
                    break;
                }

                // Transition to exhausted T cell if pd_l1 exceeds pd_1
                if (this.pd_l1 > this.pd_1) {
                    this.type = exhaustedTcells;
                    break;
                } else {
                    // Killing logic
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && G.GetAgent(movdivHood[j]).type == Myeloma && this.hours_since_kill >= 1) {
                            G.GetAgent(movdivHood[j]).Tcell_Kill();
                            this.pd_l1 += 1;
                            this.hours_since_kill = 0; // Reset kill timer

                            // Division logic
                            if (this.hours_since_division >= 24 && emptyNeighbors > 0) {
                                for (int i = 0; i < emptyNeighbors; i++) {
                                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                                    int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                        agentLattice child = G.NewAgentSQ(chosenCell);
                                        child.type = this.type;
                                        child.SetCellColor();
                                        child.tcellAge = 0;
                                        child.pd_1 = Math.floor(this.pd_1 / 4);
                                        child.pd_l1 = this.pd_l1;
                                        child.hours_since_division = 0;
                                        this.hours_since_division = 0;
                                        break;
                                    }
                                }
                            }
                            encounteredMyeloma = true; // Set the flag to true if a myeloma cell is encountered
                            break; // Exit the inner loop if a myeloma cell is encountered
                        }
                    }
                }

                // If a myeloma cell is encountered, break out of the outer loop
                if (encounteredMyeloma) {
                    break;
                }

                // Movement logic
                int moveToIndex = seekCXCL9(); // Explicitly define since seekTGFB has random values
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }


        if (type == Myeloma) {
            if (this.BCMA== false){
                this.color = RGB(100, 0, 100);

            }
            //STEP 1: death
            if (G.rn.Double() < ProbScale(G.Myeloma_DEATH_PROB, 1)) {
                Dispose();
                return;
            }
            //STEP 2: division
            if (G.rn.Double() < ProbScale(G.Myeloma_PROLIFERATION_PROB, 1)) {
                int[] divHood = MooreHood(true);
                int emptyNeighbors = MapEmptyHood(divHood);
                //Create new Agent
                if (emptyNeighbors > 0) {
                    agentLattice child = G.NewAgentSQ(divHood[G.rn.Int(emptyNeighbors)]);
                    if(G.rn.Double() < mutProb){
                        child.BCMA = false;
                    }
                    else if (G.rn.Double() > mutProb){
                        child.BCMA = true;
                    }
                    child.type = this.type;
                    child.SetCellColor();
                }
            }
        }

        if (type == inactiveTcell) {
            boolean agedDeath = false; //Tracking for aged Tcells
            if (timeStep % 24.0== 0) {
                this.tcellAge += 1;
            }
            if(G.TCE_ON){
                this.type = activeTcell;
            }
            for (int run = 0; run < 3; run++) {
                double rn_BirthDeath = G.rn.Double();
                //double pdiv = G.T_CELL_DIV_RATE;
                int[] movdivHood = MooreHood(true); // for division and movement
                int options = MapOccupiedHood(movdivHood); // mapping occupied spots
                int emptyNeighbors = MapEmptyHood(movdivHood); // mapping empty spots
                if (this.tcellAge == 30) {
                    if (emptyNeighbors > 0) {
                        for (int i = 0; i < emptyNeighbors; i++) {
                            int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                            int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                            if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                agentLattice child = G.NewAgentSQ(chosenCell);
                                child.type = this.type;
                                child.tcellAge = 0;
                                child.pd_1 = 0;
                                child.pd_l1 = 0;
                                child.hours_since_division = 0;
                                child.SetCellColor();
                                break;
                            }
                        }
                        agedDeath = true;
                        this.Dispose();
                    }
                }
                if (agedDeath){
                    break;
                }
                if (emptyNeighbors > 0) {
                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                    int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                        MoveSQ(chosenCell);
                    }
                }
            }
        }
    }
}
