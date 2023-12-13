var pointCount = 350; //  100 are 8.33 hours
var i, r;

/* - Weakly transcribed genes: 0.1 to 1 transcripts per minute (6 to 60 transcripts per hour)
- Moderately transcribed genes: 1 to 10 transcripts per minute (60 to 600 transcripts per hour)
- Highly transcribed genes: greater than 10 transcripts per minute (over 600 transcripts per hour)


- **Rapidly degraded proteins:**  Degradation rate: >1 h^-1  Half-life: <1 hour
- **Stably expressed proteins:**   Degradation rate: between 0.1 and 1 h^-1   Half-life: anywhere from 1 to 10 hours
- **Very stable proteins:**   Degradation rate: <0.1 h^-1   Half-life: >10 hours

- **mRNA Degradation rate:** on the order of 0.5 to 3 min^-1 (30 to 180 h^-1)


- **High-affinity binding:** K_d values are in the nanomolar range (nM, 10^-9 M). This indicates very tight binding, and such repressors or activators will be strongly associated with their DNA targets even at low concentrations.
- **Moderate-affinity binding:** K_d values are in the low micromolar range (µM, 10^-6 M). Proteins with moderate binding affinities may require higher concentrations to efficiently bind to their DNA targets.
- **Low-affinity binding:** K_d values are in the high micromolar to millimolar range (mM, 10^-3 M). Such interactions are weaker, and proteins may dissociate more readily from the DNA, potentially regulating target genes in a more transient fashion
 */

let prod_R3 = 3.0; // (3.33 minutes^-1 acc to Zhuo Chen) 1 (0.84 minutes^-1 acc to Shristopher A. Voigt) (0.84 acc to Jennifer S. Hallinan)
let prod_R1 = 3.0; //(3.33 minutes^-1 acc to Zhuo Chen) 1 (48 minutes^-1 acc to Shristopher A. Voigt)
let prod_L = 2.5; // (3.33 minutes^-1 acc to Zhuo Chen) 1 (6 minutes^-1 acc to Jennifer S. Hallinan)

let prod_MP1 = 1.0 ; // (1.416 transcript/minute acc to Zhuo Chen)1  (9 transcript/minute acc to Shristopher A. Voigt) (1.2 trancript/minute acc to Jennifer S. Hallinan)
let prod_MP3 = 1.0; // (1.667 transcript/minute acc to Zhuo Chen)1 (16.8 transcript/minute acc to Shristopher A. Voigt)
let prod_ML = 1.0; //   (2.083 transcript/minute acc to Zhuo Chen)1 (12 transcript/minute acc to Jennifer S. Hallinan)

let deg_R = 1.4e-2;  //  (3.33e-3 1/minute acc to Zhuo Chen) 1 (0.12 1/minute acc to Shristopher A. Voigt)
let deg_I = 1.4e-2; //  (3.33e-3 1/minute acc to Zhuo Chen) 1 (1.2 1/minute acc to Shristopher A. Voigt)
let deg_L = 1.4e-2; // (0.01 1/minute acc to Zhuo Chen) 1

//let deg_MP = 0.5; //  (0.138 1/minute acc to Zhuo Chen)1 

let for_com_RI =  1e-3; //  (5.33e-3 molecules^-1*min^-1 acc to Zhuo Chen)1 (9.06e-3 molecules^-1*min^-1 acc to Joseph A. Newman)
let for_com_RL = 1e-3; //  (5.33e-3 molecules^-1*min^-1 acc to Zhuo Chen)1 (1.61e-3 molecules^-1*min^-1 acc to Joseph A. Newman)

let n_R = 4; 
let n_A = 4; 
let K_R = 110; // (476 - 964 molecules acc Joseph A. Newman)
let K_A = 100; 



let init_L=0.9705701//664; // at the beginning, there is no SlrR. SinR is inhibiting the transcription of SlrR.
let init_I=0//24; // at the beginning, there is no SinI. AbrB is inhibiting the transcription of SinI and there's no Spo0A-P to activate it.
let init_R=200.3932; //<<<<####################################################


let data1 = [];
let R = [];
let I = [];
let L = [];
var data2 = [];
 R[0] = init_R;
 I[0] = init_I;
 L[0] = init_L;



   
    var c = [];
    var time_step = 0.11;

    for(i = 0; i < pointCount; i++) {

        if (i<150) {
        var P =  i*2 ;  // from 0 to 300
        } else if (i>=150) {
            var P =  600 - i*2 ;  //from 300 to 0
        }
        /*  
        else if (i>=1000){
            var P =  i - 1000 ;  //from 0 to 500 again
        } 
  */
        if (P<0) {
            P=0;
        }

        let activation_L = 1 / (1 + Math.pow(R[i] / K_R, n_R));
        
        let activation_P1 = activation_L * Math.pow(P / K_A, n_A) / (1 + Math.pow(P / K_A, n_A));
       
        
        // dMP3/dt equation

        let dMP3 = 1 //- deg_MP*MP3[i];
        //MP3.push(MP3[i] + dMP3/time_step);

        // dM1/dt equation
        let dMP1 = activation_P1 //- deg_MP*MP1[i];
        //MP1.push(MP1[i] + dMP1/time_step);

        // dML/dt equation
        let dML = activation_L //- deg_MP*ML[i];
        //ML.push(ML[i] + dML/time_step);

        // dR/dt equation
        let dR = ( dMP3)*prod_R3 - deg_R * R[i]- for_com_RI * R[i] * I[i] - for_com_RL * R[i] * L[i]; 
       
       
       // if (R[i] + dR/time_step < 0) {
       //     R.push(0);
      //  } else {
        R.push(R[i] + dR/time_step);
      //  }
        

        // dI/dt equation
        let dI = (dMP1)*prod_R1 - deg_I * I[i]- for_com_RI* R[i] * I[i];           
        // ; /* + dis_com_I * ComI[i]; */

      //  if (I[i] + dI/time_step < 0) {
     //       I.push(0);
      //  } else {
        I.push(I[i] + dI/time_step);
      //  }

        // dL/dt equation
        
        let dL = (dML)*prod_L - deg_L * L[i]- for_com_RL * R[i] * L[i];


      //  if (L[i] + dL/time_step < 0) {
      //      L.push(0);
     //   } else {
        L.push(L[i] + dL/time_step);
      //  }

        c.push(i)
        

    }

// Create new arrays for concentrations
var R_conc = [];
var I_conc = [];
var L_conc = [];

// After all calculations
for(i = 0; i < pointCount; i++) {
    R_conc.push(R[i]);// / 2); //From molecules to nM
    I_conc.push(I[i]);// / 2);
    L_conc.push(L[i]);// / 2);
}

    data1.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: R_conc, //SinR
        y: I_conc, //SinI
        z: L_conc, //SlrR
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

/* data2.push({
    type: 'scatter3d',
    mode: 'lines+markers',
    x: MP3, //SinR
    y: MP1, //SinI
    z: ML, //SlrR
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
}); */


/* Plotly.newPlot('myDiv2', data2, {
    scene: {
        xaxis: {
            title: 'SinR mRNA [Molecules]',
            titlefont: {
                color: 'red'
            }
        },
        yaxis: {
            title: 'SinI mRNA [Molecules]',
            titlefont: {
                color: 'blue'
            }
        },
        zaxis: {
            title: 'SlrR mRNA [Molecules]',
            titlefont: {
                color: 'rgb(255, 200, 0)'
            }
        }
    },
    showlegend: false
}); */

Plotly.newPlot('myDiv', data1, {
    scene: {
        xaxis: {
            title: 'SinR [molecules/cell]',
            titlefont: {
                color: 'red'
            }
        },
        yaxis: {
            title: 'SinI [molecules/cell]',
            titlefont: {
                color: 'blue'
            }
        } ,
        zaxis: {
            title: 'SlrR [molecules/cell]',
            titlefont: {
                color: 'rgb(180, 170, 0)'
            }
        } 
    },
    showlegend: false
});


function downloadCSV(data) {
    for (var j=0; j < data.length; j++){
        var csv = 'SinR,SinI,SlrR\n';
        data[j].x.forEach(function(item, i) {
            csv += data[j].x[i] + ',' + data[j].y[i] + ',' + data[j].z[i] + '\n';
        });

        var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
        var link = document.createElement("a");
        var url = URL.createObjectURL(blob);
        link.setAttribute("href", url);
        link.setAttribute("download", "data" + j + ".csv");
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
}

// Call the function with your data
//downloadCSV(data1);
