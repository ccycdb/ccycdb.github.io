document.addEventListener("DOMContentLoaded", function () {
  const toggleContent = document.getElementById("toggle-content");
  const secondLayer = document.querySelector(".second-layer");
  const thirdLayer = document.querySelector(".third-layer");
  const fourthLayer = document.querySelector(".fourth-layer");
  const fifthLayer = document.querySelector(".fifth-layer");

  let isCollapsed = true;

  toggleContent.addEventListener("click", function () {
    if (isCollapsed) {
      secondLayer.style.display = "block";
      isCollapsed = false;
    } else {
      secondLayer.style.display = "none";
      isCollapsed = true;
    }
  });

  secondLayer.addEventListener("click", function () {
    if (isCollapsed) {
      thirdLayer.style.display = "block";
      isCollapsed = false;
    } else {
      thirdLayer.style.display = "none";
      isCollapsed = true;
    }
  });

  thirdLayer.addEventListener("click", function () {
    if (isCollapsed) {
      fourthLayer.style.display = "block";
      isCollapsed = false;
    } else {
      fourthLayer.style.display = "none";
      isCollapsed = true;
    }
  });

  fourthLayer.addEventListener("click", function () {
    if (isCollapsed) {
      fifthLayer.style.display = "block";
      isCollapsed = false;
    } else {
      fifthLayer.style.display = "none";
      isCollapsed = true;
    }
  });
});
