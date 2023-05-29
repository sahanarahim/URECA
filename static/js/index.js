function submitGeneForm(event) {
    const titleInput = document.getElementById("gene_id");
    if (titleInput.value === "") {
        titleInput.value = titleInput.placeholder;
    }
}

function submitNameForm(event) {
    const titleInput = document.getElementById("author");
    if (titleInput.value === "") {
        titleInput.value = titleInput.placeholder;
    }
}
function submitTitleForm(event) {
    const titleInput = document.getElementById("title");
    if (titleInput.value === "") {
        titleInput.value = titleInput.placeholder;
    }
}
function submitFormWith(event, value) {
    event.preventDefault();
    // Submit the form with the clicked value
    const link = event.target;
    const form = link.closest('form');
    const input = form.querySelector('input[type="text"]');
    input.value = value;
    form.submit();
}

/*
This part of the file contains JavaScript for the help part of the landing page.
*/

// Selects the necessary HTML elements to change / attach event listeners to:
const help_information = document.querySelector('#help-text');
const buttons = document.querySelector('.button-group.hollow').querySelectorAll('.button');

// Change help information depending on which button is clicked:
const change_help_text = (text) => {
    let help_text = '';
    switch (text) {
        case 'word':
            help_text = 'Finds all entities containing the searched term.'
            break;
        case 'exact':
            help_text = 'Finds entities containing the exact term (e.g., searching for "cesa" would return "cesa expression", "cesa tracking", and so on).';
            break;
        case 'alias':
            help_text = 'Finds gene aliases (e.g., searching for "cesa1" finds terms such as "rsw1", "atcesa1", and "any1").';
            break;
        case 'substring':
            help_text = 'Finds entities containing the substring (e.g., searching for "hair" would return "root hairs", "hairy roots", and so on).';
            break;
        case 'non-alphanumeric':
            help_text = '"fer" finds entities such as "fer-tor" and "fer/ripk".';
            break;
        default:
            help_text = 'Click one of the above buttons to find out more about each search function!';
            break;
    };
    help_information.innerText = help_text;
}

// Attach event listeners to each button:
buttons.forEach(button => {
    let button_text = button.innerText.toLowerCase();
    button.setAttribute('onclick', `change_help_text('${button_text}')`);
})