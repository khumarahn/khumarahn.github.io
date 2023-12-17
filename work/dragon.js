let dragon = document.getElementById("dragon");

let disp_x = 31;
let disp_y = 17;
let disp = Array(disp_y).fill(null).map(()=>Array(disp_x).fill(' '));

function disp_text() {
  let s = '';
  for (let y = 0; y < disp_y; y++) {
    for (let x = 0; x < disp_x; x++) {
      s += disp[y][x];
    }
    s += '\n';
  }
  return s;
}

function disp_add(buf, buf_x, buf_y) {
  for (let y = 0; y < buf.length; y++) {
    for (let x = 0; x < buf[y].length; x++) {
      if (buf[y][x] != ' ' 
          && buf_x + x >= 0 && buf_x + x < disp_x 
          && buf_y + y >= 0 && buf_y + y < disp_y) {
          let c = buf[y][x];
          if (buf[y][x] == '<') {
              c = '&lt;';
          } else if (buf[y][x] == '>') {
              c = '&gt;';
          }
        disp[buf_y + y][buf_x + x] = c;
      }
    }
  }
}

function disp_del(buf, buf_x, buf_y) {
  for (let y = 0; y < buf.length; y++) {
    for (let x = 0; x < buf[y].length; x++) {
      if (buf[y][x] != ' ') {
        disp[buf_y + y][buf_x + x] = ' ';
      }
    }
  }
}

function on_resize() {
  let width_em = window.innerWidth / parseFloat(getComputedStyle(document.querySelector('body'))['font-size']);
  if (width_em > 72) {
    dragon.style.position = "fixed";
  } else {
    dragon.style.position = "static";
  }
}
on_resize();
window.addEventListener("resize", on_resize, false);

let wings_x = 0;
let wings_y = 0;
let wings = [[
'                            '.split(''),
'                            '.split(''),
'        /\\\\____             '.split(''),
'       /       \\\\           '.split(''),
'      /  / /    \\\\/         '.split(''),
'     /  / / / __            '.split(''),
'    /  / / / /              '.split(''),
'   /  / /\\/\\/               '.split(''),
'  /  \\/                     '.split(''),
'  \\//                       '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split('')
],[
'                            '.split(''),
'        /\\\\____             '.split(''),
'       /      \\\\            '.split(''),
'      /  / /   \\\\__         '.split(''),
'     /  / / / ___           '.split(''),
'    /  / / / /  /           '.split(''),
'   /  / /\\/\\/               '.split(''),
'  /  \\/                     '.split(''),
'  \\//                       '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split('')
]];

let head_x = 17;
let head_y = 0;
let head=[[
' /\\/\\___   '.split(''),
'<  õ   *\\  '.split(''),
'<   _v_v/  '.split(''),
' |  |      '.split('')
],[
' /\\/\\__    '.split(''),
'<  V  *\\   '.split(''),
'<    vv/   '.split(''),
' |  |      '.split('')
],[
'  /\\/\\_    '.split(''),
' < õ * >   '.split(''),
' < \\_^/>   '.split(''),
'  |  |     '.split('')
],[
'  /\\_/\\    '.split(''),
' <õ * õ>   '.split(''),
' <\\_^_/>   '.split(''),
'  |  |     '.split('')
],[
'  /\\_/\\    '.split(''),
' <õ * õ>   '.split(''),
' <\\_V_/>   '.split(''),
'  |  |     '.split('')
]];


let hands_x = 19;
let hands_y = 3;
let hands=[[
'         '.split(''),
'   \\     '.split(''),
'  \\_\\\\__ '.split(''),
'\\_____ \\\\'.split(''),
'      JJJ'.split(''),
],[
'     ___  '.split(''),
'   \\/  \\\\ '.split(''),
'  \\/ /JJJ '.split(''),
'\\___/     '.split(''),
'          '.split('')
]];

let body_x = 0;
let body_y = 0;
let body = [
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'                            '.split(''),
'               /            '.split(''),
'             _/     ~~~)    '.split(''),
'           _/       \\~~)    '.split(''),
'          /         /~/     '.split(''),
'             __\\   / /      '.split(''),
'         ___/  |  /_/_      '.split(''),
'               (__  \\ \\     '.split(''),
'                  JJJJJJ    '.split('')
];

let tail_x = 0;
let tail_y = 9;
let tail = [[
'         _'.split(''),
'       _/ '.split(''),
'  ____/   '.split(''),
'  \\_____/ '.split(''),
'          '.split(''),
],[
'         _'.split(''),
'   __  _/ '.split(''),
'   \\ \\/   '.split(''),
'    \\___/ '.split(''),
'          '.split(''),
],[
'         _'.split(''),
'       _/ '.split(''),
'  ____/   '.split(''),
'  \\_____/ '.split(''),
'          '.split(''),
],[
'__        '.split(''),
'\\ \\______/'.split(''),
' \\________'.split(''),
'          '.split(''),
'          '.split(''),
],[
'         _'.split(''),
'       _/ '.split(''),
'  ____/   '.split(''),
'  \\_____/ '.split(''),
'          '.split(''),
]];

let click_x = 0;
let click_y = 0;
let click = [
    ['   (CLICK)   '.split('')],
    ['  ((CLICK))  '.split('')],
    [' (( CLICK )) '.split('')]
];

// bottom text
let btext_x = 0;
let btext_y = 14;
let btext = [
    'dragon learning'.split(''),
    '       true meaning of chaos'.split('')
];

let seq = 0;

let seq_wings = 0;
let seq_wings_p = 0;

let seq_head = 0;
let seq_head_p = 0;
let seq_head_up = 0;

let seq_hands = 0;
let seq_hands_p = 0;

let seq_tail = 0;
let seq_tail_p = 0;

let seq_click = 0;
let seq_click_p = 0;

function drag() {
  if (seq == 0 || seq_wings_p != seq_wings || seq_head_p != seq_head || seq_hands_p != seq_hands || seq_tail_p != seq_tail) {
      
    disp_del(wings[seq_wings_p], wings_x, wings_y);
    disp_del(head[seq_head_p], head_x, head_y);
    disp_del(hands[seq_hands_p], hands_x, hands_y);
    disp_del(tail[seq_tail_p], tail_x, tail_y);

    // head *after* wings
    disp_add(wings[seq_wings], wings_x, wings_y);
    disp_add(head[seq_head], head_x, head_y);
    disp_add(hands[seq_hands], hands_x, hands_y);
    disp_add(tail[seq_tail], tail_x, tail_y);
    disp_add(body,body_x, body_y);

    // bottom text
    disp_add(btext, btext_x, btext_y);
  } 

  if (seq_click_p != seq_click) {
    disp_del(click[seq_click_p], click_x, click_y);
    disp_add(click[seq_click],   click_x, click_y);
  }


  seq_wings_p = seq_wings;
  seq_head_p  = seq_head;
  seq_hands_p = seq_hands;
  seq_tail_p  = seq_tail;
  seq_click_p = seq_click;

  if (seq % 8 == 0) {
      seq_click = Math.floor(click.length * Math.random());
  }

  if (seq_wings == 0) {
    if (Math.random() < 1./32.) {
      seq_wings = 1;
    }
  } else {
    if (Math.random() < 1./8.) {
      seq_wings = 0;
    }
  }
  
  if (seq_hands == 0) {
    if (Math.random() < 1./16.) {
      seq_hands = 1;
    }
  } else {
    if (Math.random() < 1./4.) {
      seq_hands = 0;
    }
  }
  
  if (seq_head_up == 0) {
    if (seq_head == 0 && Math.random() < 1./32.) {
      seq_head_up = 1;
    } else if (Math.random() < 1./16.) {
      seq_head_up = -1;
    }
  }

  if (seq_head_up == 1) {
    if (seq_head + 1 < head.length) {
      seq_head++;
    } else {
      seq_head_up = 0;
    }
  } else if (seq_head_up == -1) {
    if (seq_head > 0) {
      seq_head--;
    } else {
      seq_head_up = 0;
    }
  }

  if (seq_tail == 0) {
    if (Math.random() < 1./32.) {
      seq_tail = Math.floor(Math.random() * tail.length);
    }
  } else {
    if (Math.random() < 1./8.) {
      seq_tail = Math.floor(Math.random() * tail.length);
    }
  } 
 
  seq++;
  
  let new_disp_text = disp_text();
  if (dragon.innerHTML != new_disp_text) {
    dragon.innerHTML = disp_text();
  }
  
  setTimeout(drag, 400);
}

drag();
