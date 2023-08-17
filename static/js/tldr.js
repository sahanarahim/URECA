document.addEventListener("DOMContentLoaded", function () {
  // Your code here
  document
    .getElementById("chat-form")
    .addEventListener("submit", function (event) {
      event.preventDefault();
      let userInput = document.getElementById("user-input");
      let loadingDiv = document.getElementById("loading");
      let warningDiv = document.getElementById("warning");
      if (userInput.value === "") {
        userInput.value = "CESA1";
      }
      let url = `/tldr/search?user_input=${encodeURIComponent(
        userInput.value
      )}`;
      document.getElementById("submitBtn").setAttribute("disabled", "disabled");
      loadingDiv.style.display = "block";
      fetch(url)
        .then((response) => response.json())
        .then((data) => {
          let content = data.content;
          let warning = data.warning;
          if (warning !== "") {
            warningDiv.style.display = "block";
            warningDiv.innerHTML = "<p>" + warning + "</p>";
          } else {
            warningDiv.style.display = "none";
          }
          let resultDiv = document.getElementById("result");
          const regex = /\((\d+)\)/g;
          let matches = content.match(regex);
          if (matches) {
            for (let i = 0; i < matches.length; i++) {
              let number = matches[i].match(/\d+/)[0];
              content = content.replace(
                matches[i],
                '<a href="http://www.ncbi.nlm.nih.gov/pubmed/' +
                  number +
                  'target="_blank">' +
                  matches[i] +
                  "</a>"
              );
            }
          }
          resultDiv.innerHTML = content;
          document.getElementById("submitBtn").removeAttribute("disabled");
          loadingDiv.style.display = "none";
        })
        .catch((error) => {
          console.error("Error fetching response:", error);
          document.getElementById("submitBtn").removeAttribute("disabled");
          loadingDiv.style.display = "none";
        });
    });
});
