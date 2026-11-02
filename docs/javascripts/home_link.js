// Make the site name ("LibMAT") in the header and at the top of the sidebar a
// link to the home page, like the logo next to it.
document$.subscribe(function () {
  var logo = document.querySelector(".md-header__button.md-logo");
  if (!logo) return;
  var home = logo.getAttribute("href");

  function link(text) {
    var a = document.createElement("a");
    a.className = "libmat-home";
    a.href = home;
    a.textContent = text;
    return a;
  }

  var header = document.querySelector(".md-header__topic:first-child .md-ellipsis");
  if (header && !header.querySelector(".libmat-home")) {
    header.replaceChildren(link(header.textContent.trim()));
  }

  var title = document.querySelector(".md-nav--primary > .md-nav__title");
  if (title && !title.querySelector(".libmat-home")) {
    title.childNodes.forEach(function (node) {
      if (node.nodeType === Node.TEXT_NODE && node.textContent.trim()) {
        title.replaceChild(link(node.textContent.trim()), node);
      }
    });
  }
});
