var pointCount = 5000;
var i, r;

let prod_R3 = 0.84; // proteins/(transcript*minute)
let prod_R1 = 0.9; // proteins/(transcript*minute)
let prod_L = 0.8; // proteins/(transcript*minute)
let prod_MP1 = 9.0; // transcript/minute
let prod_MP3 = 10.8; // transcript/minute 
let prod_ML = 9.0; // transcript/minute
let deg_R = 3.33e-3;  // 1/minute  ## fixed
let deg_I = 3.33e-3; // 1/minute ## fixed
let deg_L = 0.01; // 1/minute ## fixed
let deg_MP = 0.3; // 1/minute

let for_com_RI =  8.17e-4; // 1/(minute*protein)
let for_com_RL = 0.0; 

let n_R = 2; 
let n_A = 0; 
let K_R = 100; 
let K_A = 0; 
let A = 6;


var initialValues = [280.0, 40.0, 10.0, 0.0, 0.0, 0.0];
var grid_size = 2;
var range = [0, 10]; // Adjust as needed
var data = [];



    var R = [initialValues[0]];
    var I = [initialValues[1]];
    var L = [initialValues[2]];
    var MP1 = [initialValues[3]];
    var MP3 = [initialValues[4]];
    var ML = [initialValues[5]];
    var c = [];

    for(i = 0; i < pointCount; i++) {

        let activation_P1 =0.5;// 1 / (1 + Math.pow(R[i] / K_R, n_R)) * Math.pow(A[i] / K_A, n_A) / (1 + Math.pow(A[i] / K_A, n_A));
       
        let activation_L = 1 / (1 + Math.pow(R[i] / K_R, n_R));
        // dMP3/dt equation

        let dMP3 = prod_MP3 - deg_MP*MP3[i];
        MP3.push(MP3[i] + dMP3/60);

        // dM1/dt equation
        let dMP1 = prod_MP1*activation_P1 - deg_MP*MP1[i];
        MP1.push(MP1[i] + dMP1/60);

        // dML/dt equation
        let dML = prod_ML*activation_L - deg_MP*ML[i];
        ML.push(ML[i] + dML/60);

        // dR/dt equation
         let dR = (MP1[i]+MP3[i])*prod_R3 - deg_R * R[i]- for_com_RI * R[i] * I[i];// - for_com_RL * R[i] * L[i]; /* + dis_com_L * ComL[i] + dis_com_I * ComI[i]; */
        R.push(R[i] + dR/60);

        
        // dI/dt equation
        let dI = (MP1[i])*prod_R1 - deg_I * I[i]- for_com_RI* R[i] * I[i];           
        // ; /* + dis_com_I * ComI[i]; */
        I.push(I[i] + dI/60);
        
        // dL/dt equation
        
        let dL = (ML[i])*prod_L - deg_L * L[i];
        L.push(L[i] + dL/60);

        c.push(i)
    }

    data.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: R, //SinR
        y: I, //SinI
        z: L, //SlrR
        line: {
            width: 5,
            color: c,
            colorscale: "Viridis",
            cmin: 0,
            cmax: pointCount
        },
        marker: {
            size: 1,
            color: c,
            colorscale: "Viridis",
            cmin: 0,
            cmax: pointCount
        }
    });


Plotly.newPlot('myDiv', data, {
    scene: {
        xaxis: {title: 'SinR'},
        yaxis: {title: 'SinI'},
        zaxis: {title: 'SlrR'}
    },
    showlegend: false
});


