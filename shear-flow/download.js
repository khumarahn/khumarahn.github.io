" use strict";

function downloadText(filename,text){
    // Set up the link
    let link = document.createElement("a");
    link.setAttribute("target","_blank");
    if(Blob !== undefined) {
        let blob = new Blob([text], {type: "text/plain"});
        link.setAttribute("href", URL.createObjectURL(blob));
    } else {
        link.setAttribute("href","data:text/plain," + encodeURIComponent(text));
    }
    link.setAttribute("download",filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}
