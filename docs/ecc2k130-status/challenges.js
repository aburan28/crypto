(function () {
  "use strict";

  var COLLECTION = {
    live: "Collecting on this site",
    local: "Local walk, no upload",
    catalog: "Listed, no walk here",
    none: "No store on this site"
  };

  function el(tag, text, className) {
    var node = document.createElement(tag);
    if (className) node.className = className;
    if (text != null) node.textContent = text;
    return node;
  }

  function row(challenge) {
    var tr = document.createElement("tr");
    tr.appendChild(el("th", challenge.name));
    tr.lastChild.scope = "row";
    tr.appendChild(el("td", challenge.field, "mono"));
    tr.appendChild(el("td", challenge.kind));
    var state = el("td", COLLECTION[challenge.collection] || challenge.collection);
    if (challenge.status_href) {
      var link = el("a", "Open the live counts");
      link.href = challenge.status_href;
      state.appendChild(document.createTextNode(" "));
      state.appendChild(link);
    }
    tr.appendChild(state);
    tr.appendChild(el("td", challenge.join));
    return tr;
  }

  function render(catalog, groupId) {
    var body = document.getElementById("challenge-body");
    var note = document.getElementById("group-note");
    while (body.firstChild) body.removeChild(body.firstChild);
    var shown = 0;
    catalog.groups.forEach(function (group) {
      if (groupId !== "all" && group.id !== groupId) return;
      if (note && groupId !== "all") note.textContent = group.summary;
      group.challenges.forEach(function (challenge) {
        body.appendChild(row(challenge));
        shown += 1;
      });
    });
    if (groupId === "all" && note) {
      note.textContent = catalog.claim;
    }
    if (shown === 0) {
      var empty = document.createElement("tr");
      var cell = el("td", "No challenges in this group.");
      cell.colSpan = 5;
      empty.appendChild(cell);
      body.appendChild(empty);
    }
  }

  function filters(catalog) {
    var bar = document.getElementById("challenge-filters");
    var options = [{ id: "all", title: "All" }].concat(catalog.groups.map(function (group) {
      return { id: group.id, title: group.title };
    }));
    options.forEach(function (option, index) {
      var button = el("button", option.title, "mode-tab");
      button.type = "button";
      button.setAttribute("aria-pressed", index === 0 ? "true" : "false");
      button.addEventListener("click", function () {
        Array.prototype.forEach.call(bar.children, function (child) {
          child.setAttribute("aria-pressed", "false");
        });
        button.setAttribute("aria-pressed", "true");
        render(catalog, option.id);
      });
      bar.appendChild(button);
    });
  }

  fetch("./challenges.json", { cache: "no-store" })
    .then(function (response) {
      if (!response.ok) throw new Error("challenges.json " + response.status);
      return response.json();
    })
    .then(function (catalog) {
      document.getElementById("challenge-banner").textContent = catalog.claim;
      document.getElementById("challenge-banner").className = "banner mac-current";
      filters(catalog);
      render(catalog, "all");
    })
    .catch(function () {
      var banner = document.getElementById("challenge-banner");
      banner.textContent = "The challenge roster did not load. Read challenges.json directly.";
      banner.className = "banner UNREACHABLE";
    });
}());
