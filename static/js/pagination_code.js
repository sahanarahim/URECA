const results_container = document.querySelector('#body_results');
const results = (function() {
    let data = [];
    results_container.childNodes.forEach(node => {
        let text_data = node.textContent.trim().replaceAll("\n", "");
        if (text_data.length) {
            data.push(text_data.split('  ').filter(term => term.length).map(term => term.trim()))
        } 
    })
    return data;
})();
const pagination_numbers = document.querySelector('.pagination.text-center');

// We'll just roll with 10 items per page in the pagination for now, but we could 
// also have an option for users to choose how much stuff they want to show per page:

const items_per_page = 50;                                               // Maybe make this changeable via a dropdown?
const items_to_show = Math.ceil(results.length / items_per_page);        // "mock_data" is a variable that was made in the test_output_stuff.js file!

// Add the pagination numbers when the page loads:
const add_pagination_elements = () => {
    for (let i = 1 ; i <= items_to_show ; i++) {
        let element = document.createElement('li'), number = document.createElement('span')
        number.setAttribute('aria-label', `Page ${i}`);
        number.innerText = i;
        element.setAttribute('onclick', `display_items(${i})`);
        element.id = 'pagination_number';
        element.appendChild(number);

        pagination_numbers.appendChild(element);
    }
}

// END

// Displaying the active elements by the pagination element:
const display_items = (page_number) => {
    /*
    Displays items based on which pagination number is clicked.
    */
    let lower = (page_number - 1) * items_per_page, higher = items_per_page * page_number;
    show_active_elements(page_number)

    results_container.innerHTML = '';
    results.forEach((item, index) => {
        if (index < higher && index >= lower) {
            results_container.innerHTML += `<tr>` + 
            `<td> ${item[0]} </td>` + 
            `<td> ${item[1]} </td>` + 
            `<td> ${item[2]} </td>` + 
            `<td class = 'pubmed-link pubmed-hyperlink' data-pubmed-id = '${item[3]}' data-source = '${item[0]}'
            data-typa = '${item[1]}' data-target = '${item[2]}'> ${item[3]} </td>` + 
            `</tr>`;
        }
    })
}

const show_active_elements = (p) => {
    document.querySelectorAll('#pagination_number').forEach(number => {
        number.setAttribute('class', 'disabled');       // If the user clicks another number, we want to hide other selections from them.
                                                        // Change it to a more appropriate attribute later!
        let current = number.innerText;
        if (current === p.toString()) {
            number.setAttribute('class', 'current');
        }
    })
}

// END

window.addEventListener('load', () => {
    add_pagination_elements();
    display_items(1);
});