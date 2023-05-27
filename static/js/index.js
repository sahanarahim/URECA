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