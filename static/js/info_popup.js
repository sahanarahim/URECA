const paperModal = document.querySelector('.paper-overlay');      // Modal element
const TEXT_SCALE_FACTOR = 49;                                     // This looks the best..!

const addModalContent = (paperID, source, typa, target) => {
    /*
    Fetches text from the Pubmed API to be displayed on the modal component.
    */
    let pubmedURL = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${paperID}&rettype=abstract&retmode=text`;
    
    paperModal.innerHTML = `
        <div class = "modal-content" style = "font-size: ${rescaleText()}px;"> <h5> Fetching content... </h5> </div>
    `;
    paperModal.style.display = "block";
    
    fetch(pubmedURL)
        .then(function(response) {
            if (response.ok) {
                return response.text();
            }
            throw new Error('Pubmed API call unsuccessful.');
        })
        .then(function(data) {
            let sourceWords = source.split(/[,\s&+/\\']+/), targetWords = target.split(/[,\s&+/\\']+/)
            let underlinedData = data;

            for (let i = 0 ; i < sourceWords.length ; i++) {
                let sourceRegex = new RegExp(`\\b${sourceWords[i].replace(/[-\/\\^$*+?.()|[\]{}]/g, "\\$&")}\\b`, 'gi');
                underlinedData = underlinedData.replace(sourceRegex, `<span style = "color: red;">$&</span>`);
            }
            for (let i = 0 ; i < targetWords.length ; i++) {
                let targetRegex = new RegExp(`\\b${targetWords[i].replace(/[-\/\\^$*+?.()|[\]{}]/g, "\\$&")}\\b`, 'gi');
                underlinedData = underlinedData.replace(targetRegex, `<span style = "color: red;">$&</span>`);
            }

            contents = `
            <div class = "modal-content" style = "font-size: ${rescaleText()}px;">
                <h4> ${source} ${typa} ${target} </h4>
                <br> 
                <p>
                    ${underlinedData} 
                </p>
                <br> 
                <a href = ${"https://pubmed.ncbi.nlm.nih.gov/" + paperID} target = "_blank" class = "pubmed-hyperlink"> Click here </a> to view the publication.
            </div>
            `
            paperModal.innerHTML = contents;
        })
        .catch(function(error) {
            console.log(`Error while fetching data / rendering modal: ${error.message}`);
        })
}

const rescaleText = () => {
    let height = window.innerHeight, width = window.innerWidth;
    return Math.min(width / TEXT_SCALE_FACTOR, height / TEXT_SCALE_FACTOR);
}

const closeModal = () => {
    paperModal.style.display = 'none';
    paperModal.innerHTML = '';
}

let pubmedLinks = document.querySelectorAll('span.pubmed-link');
pubmedLinks.forEach(link => {
    let paperID = link.getAttribute('data-pubmed-id'), source = link.getAttribute('data-source');
    let typa = link.getAttribute('data-typa'), target = link.getAttribute('data-target');
    link.setAttribute('onclick', `addModalContent("${paperID}", "${source}", "${typa}", "${target}")`)
});

paperModal.addEventListener('click', closeModal);
